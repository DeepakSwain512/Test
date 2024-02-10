#!/usr/bin/python3

# Curvalyzer.py
#
# This program detects curvature in the temperature dependence of amide proton
# chemical shifts. Curvature detection is treated as a model selection problem.
# Nested linear and quadratic models are considered; the statistical 
# significance of the improvement of the quadratic model (over the linear) is 
# assessed using an extra sum of squares F-test. If curvature is detected, the 
# likelihood that it is due to random errors is calculated via numerical 
# simulation based on the distribution of experimentally observed residuals.
#
# Kyle Trainor, March 2018

# ------------
# Imports
# ------------

# CSV
import pandas as pd

# numerical python
import numpy as np

# scipy
import scipy, scipy.stats as stats

# to get command-line arguments
import sys

# to work with files and directories
import os

import datetime

from itertools import islice

# manipulate zip files
import zipfile

# get pathnames matching pattern
import glob


# for emailing results
import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication

# plotting
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import pylab

from celery import Celery
import time

app = Celery('curvalyzer', backend='rpc://', broker='pyamqp://guest@localhost//')

# useful constants

# 15N/1H Ξ (reference compound: liq. NH3)
n_xi = 0.10132912

# base transmitter frequencies (Hz)
bf_h = 600130000
bf_n = 60810645.0

#######################

def plot_figs(jobDir,curves,sim_curves):
  # load temperatures from header
  full_fname = os.path.join(jobDir,"Curvalyzer.csv")
  if not(os.path.isfile(full_fname)):
    log.write("File Not Found: "+full_fname+"\nCannot plot figures.\n")
    return 1
  else:
    infile = open(full_fname)
    for idx, line in enumerate(infile):
      if idx == 3:
        templist = line.strip().split(',')[1:]
        temperatures = ','.join(templist)

  ts = [temp.strip() for temp in temperatures.split(',')]

  h_cols = [t.strip()+" K 1H-RR" for t in temperatures.split(',')]
  n_cols = [t.strip()+" K 15N-RR" for t in temperatures.split(',')]

  # initialize empty dictionaries (fill from results file)
  n_columns = {}
  h_columns = {}

  if not(os.path.isfile(full_fname)):
    log.write("File Not Found: "+full_fname+"\nCannot plot figures.\n")
    return 1
  else:
    df = pd.read_csv(full_fname, skiprows=5)

    residues = df["Residue"]
    if residues.dtype == np.int64:
      numeric_ids = True
      largest_id = residues.max()
    else:
      numeric_ids = False

    for i in range(len(ts)):
      n_columns[h_cols[i].split(" K")[0]] = df[n_cols[i]]
      h_columns[n_cols[i].split(" K")[0]] = df[h_cols[i]]


    ts_int = [float(j) for j in ts]
    for i in range(len(residues)):
      print(i)
      ts_work = []
      n = []
      h = []
      if not(np.isnan(n_columns[n_cols[0].split(" K")[0]][i])) and not(np.isnan(h_columns[h_cols[0].split(" K")[0]][i])):
        for j in range(len(ts)):
          if not(np.isnan(n_columns[n_cols[j].split(" K")[0]][i])):
            n.append(n_columns[n_cols[j].split(" K")[0]][i])
          if not(np.isnan(h_columns[h_cols[j].split(" K")[0]][i])):
            h.append(h_columns[h_cols[j].split(" K")[0]][i])
            ts_work.append(ts_int[j])

        # find best linear fit
        ph = np.polyfit(ts_work,h,1)


        # lazy way to calculate residuals (hres) 
        [f,p, hres, gq] = f_calc(ts_work,h)

        f, (p1,p2) = plt.subplots(2, sharex=True, figsize=(10,10))
        plt.suptitle("Residue "+str(residues[i]))

        # plot 1H chemical shifts vs temperature
        p1.scatter(ts_work,h, s=160, facecolors='none', edgecolors='goldenrod')

        # plot 1H chemical shift vs. temperature line of best fit
        p1.plot(ts_int,np.polyval(ph, ts_int), color="black", linewidth=2)

        # plot quadratic if curvature is detected
        if i in curves:

          if p < 0.01:
            # calculate p-value 2 (likelihood that the curvature is due to random errors)
            curve_count = 0
            for sim_curve in sim_curves:
              if abs(sim_curve) >= abs(curves[i]):
                curve_count = curve_count + 1
            pv2 = curve_count/100000

            if pv2 < 0.01:
              # plot quadratic
              pqh = np.polyfit(ts_work,h,2)
              p1.plot(ts_int,np.polyval(pqh, ts_int), color="red", linewidth=2)
              if pv2 < 0.0001:
                p1.set_title("Curvature detected: p-value 1 = "+"{0:.2g}".format(p)+"; p-value 2 < 0.0001")
              else:
                p1.set_title("Curvature detected: p-value 1 = "+"{0:.2g}".format(p)+"; p-value 2 = "+"{0:.2g}".format(pv2))


        p1.set_ylabel("$\delta$ (ppm)")

        # plot residuals vs temperature
        p2.scatter(ts_work,hres, s=160, facecolors='blue', edgecolors='blue')
        p2.set_ylabel("Residual (ppm)")
        p2.set_xlabel("Temperature (K)")

        if numeric_ids:
          chars = len(str(largest_id))
          plt.savefig(jobDir+"/residue"+str(residues[i]).zfill(chars)+".png")
        else:
          plt.savefig(jobDir+'/'+str(residues[i]).strip()+".png")

        plt.close()
  return 0

  ## end of plot_figs() function definition

#######################

def email_results(jobId,projname,email):
  jobDir = 'Curvalyzer/jobs/'+jobId
  # email the results
  emailText = "***DO NOT REPLY TO THIS EMAIL***\n"
  emailText += "The 'USERNAME@gmail.com' email address is not monitored. Send questions/comments/bug reports to 'kjtrainor@uwaterloo.ca'.\n\n"
  emailText += "Results for your job '"+projname+"' can be found in the attached 'Curvalyzer.zip' archive. The main results file, named 'Curvalyzer.csv', can be opened in a spreadsheet program such as Microsoft Excel. Also included are PNG figures generated for each peak.\n\n"
  emailText += "Your unique job ID was '"+jobId+"'. Please quote this string if you email 'kjtrainor@uwaterloo.ca' with a question or bug report.\n"

  message = MIMEMultipart()
  message["Subject"] = "Shift-T Server: Curvalyzer results for "+projname
  message["From"] = "USERNAME@gmail.com"
  message["To"] = email
  message.attach(MIMEText(emailText, "plain"))

  part = MIMEBase('application', 'octet-stream')
  part.set_payload(open(os.path.join(jobDir,'Curvalyzer.zip'), 'rb').read())
  encoders.encode_base64(part)
  part.add_header('Content-Disposition', 'attachment; filename="Curvalyzer.zip"')
  message.attach(part)

  server = smtplib.SMTP_SSL(host="smtp.gmail.com", port=465)
  server.ehlo()
  server.login("USERNAME", "PASSWORD")
  server.sendmail("USERNAME@gmail.com", email, message.as_string())
  server.close()

  ## end of email_results() function definition

#######################

def email_error(jobId,projname,email):
  jobDir = 'Curvalyzer/jobs/'+jobId

  # email the results
  emailText = "***DO NOT REPLY TO THIS EMAIL***\n"
  emailText += "The 'USERNAME@gmail.com' email address is not monitored. Send questions/comments/bug reports to 'kjtrainor@uwaterloo.ca'.\n\n"
  emailText += "An error was encountered while processing your job '"+projname+"'. The log file has been sent as an attachment to this email.\n\n"
  if jobId != 'example':
    emailText += "Your unique job ID was '"+jobId+"'. Please quote this string if you email 'kjtrainor@uwaterloo.ca' with a question or bug report.\n"

  message = MIMEMultipart()
  message["Subject"] = "Shift-T Server: Curvalyzer results for "+projname
  message["From"] = "USERNAME@gmail.com"
  message["To"] = email
  message.attach(MIMEText(emailText, "plain"))

  part = MIMEBase('application', 'octet-stream')
  part.set_payload(open(os.path.join(jobDir,'logfile.txt'), 'rb').read())
  encoders.encode_base64(part)
  part.add_header('Content-Disposition', 'attachment; filename="logfile.txt"')
  message.attach(part)

  server = smtplib.SMTP_SSL(host="smtp.gmail.com", port=465)
  server.ehlo()
  server.login("USERNAME", "PASSWORD")
  server.sendmail("USERNAME@gmail.com", email, message.as_string())
  server.close()

  ## end of email_error() function definition 

#######################

def f_calc(t,cs):
  lin_fit = np.polyfit(t,cs,1)
  hres = np.polyval(lin_fit, t)-cs
  rss_lin = np.sum((hres)**2)

  quad_fit = np.polyfit(t,cs,2)
  rss_quad = np.sum((np.polyval(quad_fit, t)-cs)**2)

  f = (rss_lin-rss_quad)/(rss_quad/(len(cs)-3))
  p = 1.0 - scipy.stats.f.cdf(f,1,len(t)-2-1)
  return [f,p, hres, quad_fit[0]]

  ## end of f_calc() function definition 

#######################

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

  ## end of isFloat() function definition 

########################################
# start of function executed by Celery #
########################################

@app.task
def pyCurvalyzer(**kwargs):
  email = kwargs['email']
  jobName = kwargs['job']
  jobId = kwargs['jobId']

  # initialize empty dictionaries (fill from results file)
  n_columns = {}
  h_columns = {}


  jobDir = 'Curvalyzer/jobs/'+jobId
  jobfile = open(os.path.join(jobDir,'job.txt'), 'w')
  jobfile.write(email+'\n')
  jobfile.write(jobName+'\n')
  jobfile.write(jobId+'\n')
  jobfile.write('{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())+'\n')
  jobfile.close()

  # open log
  global log
  log = open(os.path.join(jobDir,'logfile.txt'), 'w')

  # clean up input CSV
  startDir=os.getcwd()
  os.chdir(jobDir)
  os.system('mac2unix curvalyzer_input.csv')
  os.system('dos2unix curvalyzer_input.csv')
  os.system("perl -i.bak -pe 's/[^[:ascii:]]//g' curvalyzer_input.csv")
  os.chdir(startDir)

  header=open(os.path.join(jobDir,"curvalyzer_input.csv"), 'rU')
  header_lines = list(islice(header, 4))
  projname = jobName
  references = ','.join(list(filter(None, header_lines[1].strip().split(',')[1:])))
  nom_temperatures = ','.join(list(filter(None, header_lines[2].strip().split(',')[1:])))
  temps = ','.join(list(filter(None, header_lines[3].strip().split(',')[1:])))

  # make sure that number of reference peaks matches number of temperatures
  print(references)
  print(temps)
  if len(references.split(",")) != len(temps.split(",")):
    log.write("Error: "+"The number of reference peaks must match the number of peak lists\n")
    log.close()
    email_error(jobId, jobName, email)
    return

  # make sure that the reference peaks are all numbers
  for ref in references.split(","):
    ref = ref.strip()
    if not(isFloat(ref)):
      log.write("Error: "+"One of the reference peaks is not a number\n")
      log.close()
      email_error(jobId, jobName, email)
      return
  refs = references.split(',')

  # make sure that the temperatures are all numbers
  for idx, temp in enumerate(nom_temperatures.split(",")):
    temp = temp.strip()
    if not(isFloat(temp)):
      log.write("Error: "+"One of the temperatures is not a number\n")
      log.close()
      email_error(jobId, jobName, email)
      return

  templist = list()
  for idx, temp in enumerate(temps.split(",")):
    temp = temp.strip()
    if not(isFloat(temp)):
      log.write("Error: "+"One of the temperatures is not a number\n")
      log.close()
      email_error(jobId, jobName, email)
      return
    else:
      templist.append(format(float(temp), '.1f'))

  temperatures = ', '.join(templist)

  h_cols = [t.strip()+" K 1H" for t in temperatures.split(',')]
  n_cols = [t.strip()+" K 15N" for t in temperatures.split(',')]

    
  # read chemical shift data
  df = pd.read_csv(os.path.join(jobDir,"curvalyzer_input.csv"), skiprows=5)

  # determine residue naming convention
  residues = df["Residue"]
  if residues.dtype == np.int64:
    numeric_ids = True
    largest_id = residues.max()
  else:
    numeric_ids = False

  # load chemical shifts
  for i in range(len(temperatures.split(','))):
    n_columns[h_cols[i].split(" K")[0]] = df[n_cols[i]]
    h_columns[n_cols[i].split(" K")[0]] = df[h_cols[i]]

  peaks=list(zip(df["Ass. 15N"],df["Ass. 1H"],df["Residue"]))

  # open output file handle
  csvfile = os.path.join(jobDir, 'Curvalyzer.csv')
  outhandle = open(csvfile, "w")

  # write headers
  outlist = []
  outlist.append("Project: "+projname)
  outhandle.write(",".join(outlist)+"\n")
  outlist = []
  outlist.append("Reference Peaks (ppm):")
  outlist.append(references.strip())
  outhandle.write(",".join(outlist)+"\n")
  outlist = []
  outlist.append("Nominal Temperatures (K):")
  outlist.append(nom_temperatures.strip())
  outhandle.write(",".join(outlist)+"\n")
  outlist = []
  outlist.append("Calculated Temperatures (K):")
  outlist.append(temperatures.strip())
  outhandle.write(",".join(outlist)+"\n")
  outhandle.write("\n")
  outlist = []
  outlist.append("Residue")
  outlist.append("Ass. 1H")
  outlist.append("Ass. 15N")
  outlist.append("Notes")
  outlist.append("1H Δδ/ΔT (ppb/K)")
  outlist.append("1H RSS")
  outlist.append("15N Δδ/ΔT (ppb/K)")
  outlist.append("15N RSS")
  outlist.append("") # spacer
  for x,tempx in enumerate(temperatures.split(",")):
    outlist.append(tempx.strip()+" K 1H")
    outlist.append(tempx.strip()+" K 15N")
  outlist.append("") # spacer
  for x,tempx in enumerate(temperatures.split(",")):
    outlist.append(tempx.strip()+" K 1H-RR")
    outlist.append(tempx.strip()+" K 15N-RR")
  outhandle.write(",".join(outlist)+"\n")


  ccount = 0
  tcount = 0
  scount = 0
  cres = []
  ncres = []
  curves = {}
  n_dict = {}
  h_dict = {}
  n_rr_dict = {}
  h_rr_dict = {}
  pn_dict = {}
  ph_dict = {}
  ts_work_dict = {}


  t = [float(ts) for ts in temperatures.split(",")]
  for i in range(len(residues)):
    n = []
    h = []
    h_rr = []
    n_rr = []
    ts_work = []
    #if not(np.isnan(n_columns[n_cols[0].split(" K")[0]][i])) and not(np.isnan(h_columns[h_cols[0].split(" K")[0]][i])):
    if True:
      for j in range(len(temperatures.split(','))):
        if not(np.isnan(n_columns[n_cols[j].split(" K")[0]][i])):
          n.append(n_columns[h_cols[j].split(" K")[0]][i])
        else:
          log.write("Warning: residue "+str(i+1)+" chemical shift data contains a non-numeric value.\n")
        if not(np.isnan(h_columns[h_cols[j].split(" K")[0]][i])):
          h.append(h_columns[n_cols[j].split(" K")[0]][i])
          ts_work.append(t[j])
        else:
          log.write("Warning: residue "+str(i+1)+" chemical shift data contains a non-numeric value.\n")
    if len(h) != len(n):
      log.write("Error: "+"Mismatch in 1H and 15N chemical shift data.\n")
      log.close()
      email_error(jobId, jobName, email)
      return

    if len(h) > 0:
      # recalc temperature coefficients etc.

      # rereference 
      refs = references.split(",")
      h_rr = []
      for l, h_raw in enumerate(h):
        h_rr.append(h_raw-float(refs[l]))

      n_rr = []
      for l, n_raw in enumerate(n):
        # calculate DSS frequency in Hz
        dss_freq = (float(refs[l])*(bf_h/1000000)) + bf_h
        # multiply by Ξ ratio
        n_zero_freq = dss_freq * n_xi
        # find difference between calculated zero ppm and transmitter freq. 
        n_ppm_adj = (n_zero_freq - bf_n)/bf_n*1000000
        # rereference
        n_rr.append(n_raw-n_ppm_adj)

      # calculate temperature coefficient using rereferenced shifts
      ph = np.polyfit(t[0:len(h_rr)],h_rr,1)
      pn = np.polyfit(t[0:len(n_rr)],n_rr,1)
      #phn = np.polyfit(h_rr,n_rr,1)

      # debug output
      print("RES:",peaks[i][2])

      
      # ignore short lines; otherwise, calculate residuals and sort into curved ('cres') and non ('ncres')
      if len(t)>len(h_rr):
        print("\tshort.")
        scount = scount + 1
      else:
        [f,p, hres, quad] = f_calc(t,h_rr)
        flag = 0
        for k in range(len(t)):
          [fs,ps, hres_s, quad_s] = f_calc(t[0:k]+t[k+1:],h_rr[0:k]+h_rr[k+1:])
          if ps >= 0.01:
            flag=1

        if p<0.01 and not(flag):
          print("\tcurved.")
          ccount=ccount+1
          curves[i] = quad
          cres = np.concatenate((cres,hres),axis=0)
        else:
          ncres = np.concatenate((ncres,hres),axis=0)

        tcount=tcount+1

        print()

      pn_dict[i] = pn
      ph_dict[i] = ph
      ts_work_dict[i] = ts_work
      
    n_dict[i] = n
    h_dict[i] = h
    n_rr_dict[i] = n_rr
    h_rr_dict[i] = h_rr

  # fit residuals (non-curved only!) to t-distribution
  param = stats.t.fit(ncres)
  # debug output
  print("non-curved residuals: t-distribution fit parameters:")
  print("\tdf = "+str(list(param)[0]))
  print("\tloc = "+str(param[1]))
  print("\tscale = "+str(param[2]))
  x = np.linspace(ncres.min(), ncres.max(), 100)
  
  # plot fitted t-distribution overlaid on histogram of non-curved residuals
  p = stats.t.pdf(x, loc=param[1], scale=param[2], df=param[0])
  histo, bin_edges = np.histogram(ncres, bins='auto', normed=False)
  number_of_bins = len(bin_edges)-1
  scaling_factor = len(ncres)*(ncres.max()-ncres.min())/number_of_bins
  plt.close()
  plt.plot(x, scaling_factor*p, 'k', linewidth=2)
  plt.hist(ncres, bins='auto', normed=False)
  title = "Non-curved Residuals: Histogram and t-Distribution Fit"
  plt.title(title)
  plt.savefig(jobDir+"/residuals_t-distribution.png")

  # perform numerical simulation
  #  - calculate 100000 sets of 'fake' residuals
  #    (i.e., drawn from above t-distribution)
  #  - for any sets of 'fake' residuals that test
  #    positive for curvature, record the quadratic
  #    coefficient (for comparison with those from
  #    'real' sets of residuals in the plot_figs 
  #    function)
  sim_curves = []
  for i in range(100000):

    # monitor progress (debug purposes only)
    if i%1000 == 0:
      print(i)

    # simulate residuals
    r = list(stats.t.rvs(param[0],param[1],param[2],len(temperatures.split(','))))

    # check simulated residuals for curvature
    [f, p, hres, quad] = f_calc(t,r)

    flag = 0
    for k in range(len(t)):
      [fs,ps, hres_s,q] = f_calc(t[0:k]+t[k+1:],r[0:k]+r[k+1:])
      if ps >= 0.01:
        flag=1

    if flag and p<0.01:
      sim_curves.append(quad)

  # debug output
  print(sim_curves)
  print(len(sim_curves), len(sim_curves)/100000)


  # generate per-residue CSV output and PNG figures
  for i in range(len(residues)):
    print(i)
    n = n_dict[i]
    h = h_dict[i]
    n_rr = n_rr_dict[i]
    h_rr = h_rr_dict[i]
    if len(h) > 0:
      pn = pn_dict[i]
      ph = ph_dict[i]
      ts_work = ts_work_dict[i]
      # construct a line of output in list form
      outlist = []
      outlist.append(str(peaks[i][2])) # residue identifier
      outlist.append(str(peaks[i][1])) # 1H
      outlist.append(str(peaks[i][0])) # 15N

      if len(h_rr) == len(t):
        if i in curves:
          outlist.append("Curvature detected.") # notes 
        else:
          outlist.append("") # notes (blank)
      else:
        outlist.append("Short line.") # notes 

      outlist.append(str(ph[0]*1000)) # 1H temp coefficient in ppb/K
      outlist.append(str(np.sum((np.polyval(ph, t[0:len(h_rr)]) - h_rr) ** 2))) # RSS
      outlist.append(str(pn[0]*1000)) # 15N temp coefficient in ppb/K
      outlist.append(str(np.sum((np.polyval(pn, t[0:len(n_rr)]) - n_rr) ** 2))) # RSS
      #outlist.append(str(phn[0])) # slope in the 1H - 15N plane
      #outlist.append(str(np.sum((np.polyval(phn, h_rr) - n_rr) ** 2))) # RSS

      outlist.append("") # column spacer

      for l in range(len(h)):
        outlist.append(str(h[l]))
        outlist.append(str(n[l]))

      # pad short lines so that columns in the csv output line up
      for l in range(len(temperatures.split(','))-len(h)):
        outlist.append("") # column spacer
        outlist.append("") # column spacer

      outlist.append("") # column spacer

      for l in range(len(h_rr)):
        outlist.append(str(h_rr[l]))
        outlist.append(str(n_rr[l])) 

      # pad short lines so that columns in the csv output line up
      for l in range(len(temperatures.split(','))-len(h_rr)):
        outlist.append("") # column spacer
        outlist.append("") # column spacer


    elif np.isnan(peaks[i][1]) or np.isnan(peaks[i][0]):
      outlist = []
      outlist.append(str(peaks[i][2])) # residue identifier

    else:
      outlist = []
      outlist.append(str(peaks[i][2])) # residue identifier
      outlist.append(str(peaks[i][1])) # 15N 
      outlist.append(str(peaks[i][0])) # 1H
      outlist.append("No solution found.") # notes

    outhandle.write(",".join(outlist)+"\n")


  outhandle.close()

  # debug output
  print(ccount,tcount, scount)

  if (plot_figs(jobDir,curves,sim_curves)!=0):
    log.close()
    email_error(jobId, jobName, email)
    return

  zfname = os.path.join(jobDir,'Curvalyzer.zip')
  zfile = zipfile.ZipFile(zfname,'w')
  zfile.write(csvfile, 'Curvalyzer.csv', zipfile.ZIP_DEFLATED)

  for fname in glob.glob(jobDir+"/*.png"):
    zfile.write(fname, os.path.basename(fname), zipfile.ZIP_DEFLATED)

  zfile.close()
  email_results(jobId, jobName, email)
