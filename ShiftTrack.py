#!/usr/bin/python3 

# ShiftTrack.py
#
# This is a program for the automatic calculation of temperature coefficients 
# from 1H-15N correlation spectra. It operates on the assumption that the 
# peaks corresponding to a particular amide resonance at different temperatures
# fall (approximately) on a line in the 1H-15N plane, and that there is an upper
# limit on the distance between consecutive points. A weak constraint on the
# standard deviation of the distances between the points has also been found to 
# be useful. 
#
# Kyle Trainor, March 2018

# ------------
# Imports
# ------------

# useful stuff for parsing CSV
import pandas as pd

# useful stuff for linear regression
import numpy as np

# to get command-line arguments
import sys

# to work with files and directories
import os

# to manipulate iterators
from itertools import chain

# plotting
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

# scipy
import scipy, scipy.stats as stats

from celery import Celery
import time

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

app = Celery('shifttrack', backend='rpc://', broker='pyamqp://guest@localhost//')

# ------------
# Constants
# ------------

# relevant column headers in peak data exported from TopSpin
f1 = "ν(F1) [ppm]"
f2 = "ν(F2) [ppm]"

# alternate column headers specified in FAQ
f1_alt = "15N"
f2_alt = "1H"

# ------------
# Empirically tuned parameters 
# ------------

# RSS cut-off for linearity; for best results keep this
# relaxed (i.e. not too small).
rss_cutoff = 0.5

# std dev cut-offs for point spacing in each dimension
stdev_h_cutoff = 0.015
stdev_n_cutoff = 0.100

# largest spacing outlier allowed (in std devs; 5 is very generous)
outlier_h_cutoff = 5
outlier_n_cutoff = 5


# ------------
# The heavy lifting
# ------------

def increasing(L):
    return all(x<=y for x, y in zip(L, L[1:]))

def decreasing(L):
    return all(x>=y for x, y in zip(L, L[1:]))

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

# 'main_program' called from section below
def main_program():
  sol_count = 0
  ass_count = 0
  
  # initialize empty dictionary (will contain temperature data)
  temps = {}

  # load variable temperature data into 'temps' dictionary
  files = [f.strip() for f in peaklists.split(",")]
  for i,fname in enumerate(files):
    full_fname = directory + "/" + fname
    if not(os.path.isfile(full_fname)):
      log.write("File Not Found: "+full_fname+'\n')
      return 1
      
    df = pd.read_csv(full_fname)

    if (f1 in df.axes[1]) and (f2 in df.axes[1]):
      f1ppm = df[f1]
      f2ppm = df[f2]
    elif ('15N' in df.axes[1]) and ('1H' in df.axes[1]):
      f1ppm = df['15N']
      f2ppm = df['1H']
    else:
      log.write("Error: "+full_fname+" peak list format error\n")
      return 1

    temps[fname] = list(zip(f1ppm,f2ppm))

  if not(os.path.isfile(assignments.strip())):
    log.write("File Not Found: "+assignments+'\n')
    return 1

  peaks = []
  infile = open(assignments.strip())
  for line in infile:
    splitline = line.split(",")
    # import peak list; name of assignment stored in last position of tuple
    peaks.append([splitline[1].strip(),splitline[2].strip(),splitline[0].strip()])

  # make sure that number of reference peaks matches number of temperatures
  if len(references.split(",")) != len(peaklists.split(",")):
    log.write("Error: "+"The number of reference peaks must match the number of peak lists\n")
    return 1

  # make sure that number of reference peaks matches number of temperatures
  if len(nom_temperatures.split(",")) != len(peaklists.split(",")):
    log.write("Error: "+"The number of temperatures must match the number of peak lists\n")
    return 1

  # make sure that the reference peaks are all numbers
  for ref in references.split(","):
    ref = ref.strip()
    if not(isFloat(ref)):
      log.write("Error: "+"One of the reference peaks is not a number\n")
      return 1

  # make sure that the nominal temperatures are all numbers and calculate temperatures based upon DSS peaks
  templist = []
  refs = references.split(',')
  for idx, temp in enumerate(nom_temperatures.split(",")):
    temp = temp.strip()
    if not(isFloat(temp)):
      log.write("Error: "+"One of the temperatures is not a number\n")
      return 1
    else:
      if len(templist) == 0 or calc_temp != 1:
        last_temp = float(temp)
        templist.append(format(last_temp, '.1f'))
      else:
        delta_ref = float(refs[idx])-float(refs[idx-1])
        deltaT = delta_ref/0.0119 # based upon known 11.9 ppb/deg temperature dependence of water
        last_temp = last_temp + deltaT
        templist.append(format(last_temp, '.1f'))

    temperatures = ', '.join(templist)

  # make sure that starting temperature index is a valid number
  if not(start_index.isnumeric()):
    log.write("Error: "+"Starting temperature index must be an integer\n")
    return 1
  if (int(start_index) < 1 or int(start_index) > len(temperatures.split(","))):
    log.write("Error: "+"Starting temperature index out of range\n")
    return 1

  # open output file handle
  outhandle = open(csvfile, "w")

  # write headers
  outlist = []
  outlist.append("Project: "+projname)
  outhandle.write(",".join(outlist)+"\n")
  outlist = []
  outlist.append("Reference Peaks (ppm):")
  outlist.append(references)
  outhandle.write(",".join(outlist)+"\n")
  outlist = []
  outlist.append("Nominal Temperatures (K):")
  outlist.append(nom_temperatures)
  outhandle.write(",".join(outlist)+"\n")
  outlist = []
  outlist.append("Calculated Temperatures (K):")
  outlist.append(temperatures)
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
  outlist.append("ΔδN/ΔδH")
  outlist.append("RSS")
  outlist.append("") # spacer
  for x,tempx in enumerate(temperatures.split(",")):
    outlist.append(tempx.strip()+" K 1H")
    outlist.append(tempx.strip()+" K 15N")
  outlist.append("") # spacer
  for x,tempx in enumerate(temperatures.split(",")):
    outlist.append(tempx.strip()+" K 1H-RR")
    outlist.append(tempx.strip()+" K 15N-RR")
  outhandle.write(",".join(outlist)+"\n")

  # get index of temperature/peak list to start with
  start = int(start_index)-1


  allres = [] # all residuals

  # for each assigned peak
  for peak in peaks:

    log.write("Processing assignment "+peak[2]+"...")

    if peak[0] == "" or peak[1] == "":
      log.write("no assignment provided."+'\n')
      outlist = []
      outlist.append(peak[2]) # residue identifier
      outlist.append(str(peak[1])) # 15N 
      outlist.append(str(peak[0])) # 1H
      outlist.append("") # notes (blank)
      outhandle.write(",".join(outlist)+"\n")

    else:
      log.write('\n')
      ass_count = ass_count + 1

      # find point in data from starting temperature that most closely matches 
      # the assigned peak
      min_dist = 9999.99
      min_point = 9999
      for i, point in enumerate(temps[files[start]]):
        # dist is actually distance squared; caculating the square root gains 
        # us nothing here
        dist = (point[0]-float(peak[0]))**2 + (point[1]-float(peak[1]))**2
        if dist < min_dist:
          min_dist = dist
          min_point = i

      curr = temps[files[start]][min_point]
      
      # starting point (not really a line yet)
      lines = []
      lines.append([curr])


      shortlines = []

      # process remaining temperatures:
      #   Sarting point may be in the middle; go down from there first (prepending 
      #   points to candidate lines, then up afterwards (appending points).
      for k in chain(reversed(range(0,start)), range(start+1,len(files))):
        log.write("...working on "+files[k].split(".")[0]+'\n')

        oldlines = lines
        lines = []

        # find candidate lines: previous lines extended by a point 
        # (within distance cut-off) from spectrum at the new temp
        for j, line in enumerate(oldlines):
          # start from either the first or last point in line
          if k > start: 
            curr = line[len(line)-1]
          else:
            curr = line[0]
            
          # calculate distances
          for i, point in enumerate(temps[files[k]]):
            # dist is actually distance squared; caculating the square root 
            # gains us nothing here
            dist = (point[0]-curr[0])**2 + (point[1]-curr[1])**2
            # if within cut-off, add to list
            if dist < dist_cutoff:
              if k > start:
                lines.append(line+[point])
              else:
                lines.append([point]+line)

        oldlines = lines
        lines = []

        # assess linearity of candidate lines, discard obviously nonlinear
        #  (linear fitting after processing each new temp is somewhat
        #   inefficient, but probably preferable to allowing the number 
        #   of lines under consideration to blow up)
        for line in oldlines:
          # unpack x and y coordinates
          n,h=zip(*line)

          # calculate coeff. of determination
          r2 = (np.corrcoef(h, n)[0,1])**2

          # linear regression
          p = np.polyfit(h,n,1)

          # calculate RSS 
          rss = np.sum((np.polyval(p, h) - n) ** 2)

          # keep if linear approximation is good enough
          if rss < rss_cutoff and (increasing(h) or decreasing(h)):
            lines.append(line)

        # save shorter lines for possible consideration later on
        if k < len(files)-1:
          shortlines.append(lines)

      # check out results
      i_best = 999
      rss_min = 9999

      log.write("considering "+str(len(lines))+" full lines"+'\n')
      for i,line in enumerate(lines):
        # for std dev calculation
        ndiff = []
        hdiff = []
        for j in range(1, len(line)):
          ndiff.append(line[j][0]-line[j-1][0]) 
          hdiff.append(line[j][1]-line[j-1][1]) 

        # calculate deviations from mean spacings
        hstdev = np.std(hdiff)
        hmean = np.mean(hdiff)
        hdiff_dm = [] # list of absolute differences from the mean
        for x, diff in enumerate(hdiff):
          hdiff_dm.append(abs(diff-hmean))

        nstdev = np.std(ndiff)
        nmean = np.mean(ndiff)
        ndiff_dm = [] # list of absolute differences from the mean
        for x, diff in enumerate(ndiff):
          ndiff_dm.append(abs(diff-nmean))

        # unpack x and y coordinates
        n,h=zip(*line)
        r2 = (np.corrcoef(h, n)[0,1])**2
        p = np.polyfit(h,n,1)
        rss = np.sum((np.polyval(p, h) - n) ** 2)
        
        # if this set of points is more linear than the previous best,
        # and isn't weeded out by standard deviation cut-offs or 
        # outlier checks, it becomes the leading candidate
        if  rss < rss_min and \
            hstdev < stdev_h_cutoff and \
            nstdev < stdev_n_cutoff and \
            max(hdiff_dm) <= outlier_h_cutoff * hstdev and \
            max(ndiff_dm) <= outlier_n_cutoff * nstdev:
          best_i = i
          rss_min = rss
          r2_best = r2
          std_best_h = hstdev
          std_best_n = nstdev
          hmaxdm = max(hdiff_dm)
          nmaxdm = max(ndiff_dm)

      # if no full length lines (i.e. a point at each temperature) 
      # were found, consider shorter solutions 
      k = len(shortlines)-1
      while rss_min == 9999 and k >= (short_line_limit-1):
        lines = shortlines[k]

        log.write("considering "+str(len(lines))+" of length "+str(k+2)+'\n')
        for i,line in enumerate(lines):
          # for std dev calculation
          ndiff = []
          hdiff = []
          for j in range(1, len(line)):
            ndiff.append(line[j][0]-line[j-1][0]) 
            hdiff.append(line[j][1]-line[j-1][1]) 

          # calculate deviations from mean spacings
          hstdev = np.std(hdiff)
          hmean = np.mean(hdiff)
          hdiff_dm = [] # list of absolute differences from the mean
          for x, diff in enumerate(hdiff):
            hdiff_dm.append(abs(diff-hmean))

          nstdev = np.std(ndiff)
          nmean = np.mean(ndiff)
          ndiff_dm = [] # list of absolute differences from the mean
          for x, diff in enumerate(ndiff):
            ndiff_dm.append(abs(diff-nmean))

          # unpack x and y coordinates
          n,h=zip(*line)
          r2 = (np.corrcoef(h, n)[0,1])**2
          p = np.polyfit(h,n,1)
          rss = np.sum((np.polyval(p, h) - n) ** 2)
          
          # if this set of points is more linear than the previous best,
          # and isn't weeded out by standard deviation cut-offs or 
          # outlier checks, it becomes the leading candidate
          if  rss < rss_min and \
              hstdev < stdev_h_cutoff and \
              nstdev < stdev_n_cutoff and \
              max(hdiff_dm) <= outlier_h_cutoff * hstdev and \
              max(ndiff_dm) <= outlier_n_cutoff * nstdev:
            best_i = i
            rss_min = rss
            r2_best = r2
            std_best_h = hstdev
            std_best_n = nstdev
            hmaxdm = max(hdiff_dm)
            nmaxdm = max(ndiff_dm)

        # decrement k to look at even shorter lines if the while loop doesn't terminate
        k = k - 1



      # if a good line was found, process and write output file
      if rss_min != 9999:
        log.write("done."+'\n')
        log.flush()
        # unpack 1H and 15N coordinates
        n,h=zip(*lines[best_i])

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
        t = [float(ts) for ts in temperatures.split(",")]
        ph = np.polyfit(t[0:len(h_rr)],h_rr,1)
        pn = np.polyfit(t[0:len(n_rr)],n_rr,1)
        phn = np.polyfit(h_rr,n_rr,1)

        # calculate residuals for the purpose of estimating the standard deviation 
        # (used later to help flag points for review)
        allres = np.concatenate((allres, (np.polyval(ph,t[0:len(h_rr)])-h_rr)), axis=0)
        
        # construct a line of output in list form
        outlist = []
        outlist.append(peak[2]) # residue identifier
        outlist.append(str(peak[1])) # 1H
        outlist.append(str(peak[0])) # 15N

        if len(h_rr) == len(t):
          outlist.append("") # notes (blank)
        else:
          outlist.append("Short line.") # notes 

        outlist.append(str(ph[0]*1000)) # 1H temp coefficient in ppb/K
        outlist.append(str(np.sum((np.polyval(ph, t[0:len(h_rr)]) - h_rr) ** 2))) # RSS
        outlist.append(str(pn[0]*1000)) # 15N temp coefficient in ppb/K
        outlist.append(str(np.sum((np.polyval(pn, t[0:len(n_rr)]) - n_rr) ** 2))) # RSS
        outlist.append(str(phn[0])) # slope in the 1H - 15N plane
        outlist.append(str(np.sum((np.polyval(phn, h_rr) - n_rr) ** 2))) # RSS

        outlist.append("") # column spacer

        for l in range(len(h)):
          outlist.append(str(h[l]))
          outlist.append(str(n[l]))

        # pad short lines so that columns in the csv output line up
        for l in range(len(files)-len(h)):
          outlist.append("") # column spacer
          outlist.append("") # column spacer

        outlist.append("") # column spacer

        for l in range(len(h_rr)):
          outlist.append(str(h_rr[l]))
          outlist.append(str(n_rr[l])) 

        # pad short lines so that columns in the csv output line up
        for l in range(len(files)-len(h_rr)):
          outlist.append("") # column spacer
          outlist.append("") # column spacer


        # output to file
        outhandle.write(",".join(outlist)+"\n")

        sol_count = sol_count + 1

      else:
        log.write("no solution found."+'\n')
        outlist = []
        outlist.append(peak[2]) # residue identifier
        outlist.append(str(peak[1])) # 15N 
        outlist.append(str(peak[0])) # 1H
        outlist.append("No solution found.") # notes 
        outhandle.write(",".join(outlist)+"\n")

      log.write('\n')


  log.write('\n')
  log.write("Residues: "+str(len(peaks))+'\n')
  log.write("Assignments: "+str(ass_count)+'\n')
  log.write("Solved: "+str(sol_count)+'\n')
  outhandle.close()

  mu, stddev = stats.norm.fit(allres)
  if (plot_figs(stddev)==1):
    return 1

  zfile = zipfile.ZipFile('ShiftTrack.zip','w')
  zfile.write('logfile.txt', 'logfile.txt', zipfile.ZIP_DEFLATED)

  for fname in glob.glob("*.csv"):
    zfile.write(fname, os.path.basename(fname), zipfile.ZIP_DEFLATED)

  for fname in glob.glob("*.png"):
    zfile.write(fname, os.path.basename(fname), zipfile.ZIP_DEFLATED)

  zfile.close()

  return 0

  ## end of main_program() function definition  

#######################

def email_results():
  # email the results
  emailText = "***DO NOT REPLY TO THIS EMAIL***\n"
  emailText += "The 'USERNAME@gmail.com' email address is not monitored. Send questions/comments/bug reports to 'kjtrainor@uwaterloo.ca'.\n\n"
  emailText += "Results for your job '"+projname+"' can be found in the attached 'ShiftTrack.zip' archive. The main results file, named 'ShiftTrack.csv', can be opened in a spreadsheet program such as Microsoft Excel. Also included are PNG figures generated for each peak that was successfully followed over temperature.\n\n"
  if jobId != 'example':
    emailText += "Your unique job ID was '"+jobId+"'. Please quote this string if you email 'kjtrainor@uwaterloo.ca' with a question or bug report.\n"

  message = MIMEMultipart()
  message["Subject"] = "Shift-T Server: ShiftTrack results for "+projname
  message["From"] = "USERNAME@gmail.com"
  message["To"] = email
  message.attach(MIMEText(emailText, "plain"))

  part = MIMEBase('application', 'octet-stream')
  part.set_payload(open('ShiftTrack.zip', 'rb').read())
  encoders.encode_base64(part)
  part.add_header('Content-Disposition', 'attachment; filename="ShiftTrack.zip"')
  message.attach(part)

  server = smtplib.SMTP_SSL(host="smtp.gmail.com", port=465)
  server.ehlo()
  server.login("USERNAME", "PASSWORD")
  server.sendmail("USERNAME@gmail.com", email, message.as_string())
  server.close()

  ## end of email_results() function definition  

#######################

def email_error():
  # email the results
  emailText = "***DO NOT REPLY TO THIS EMAIL***\n"
  emailText += "The 'USERNAME@gmail.com' email address is not monitored. Send questions/comments/bug reports to 'kjtrainor@uwaterloo.ca'.\n\n"
  emailText += "An error was encountered while processing your job '"+projname+"'. The log file has been sent as an attachment to this email.\n\n"
  if jobId != 'example':
    emailText += "Your unique job ID was '"+jobId+"'. Please quote this string if you email 'kjtrainor@uwaterloo.ca' with a question or bug report.\n"

  message = MIMEMultipart()
  message["Subject"] = "Shift-T Server: ShiftTrack results for "+projname
  message["From"] = "USERNAME@gmail.com"
  message["To"] = email
  message.attach(MIMEText(emailText, "plain"))

  part = MIMEBase('application', 'octet-stream')
  part.set_payload(open('logfile.txt', 'rb').read())
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

def plot_figs(stddev):
  # load temperatures from header
  full_fname = directory + "/" + csvfile
  if not(os.path.isfile(full_fname)):
    log.write("File Not Found: "+full_fname+"\nCannot plot figures.\n")
    return 1
  else:  
    infile = open(full_fname)
    for idx, line in enumerate(infile):
      if idx == 3:  
        templist = line.strip().split(',')[1:]
        temperatures = ','.join(templist)

  files = [f.strip() for f in peaklists.split(',')]

  refs = [ref.strip() for ref in references.split(',')]
  ts = [temp.strip() for temp in temperatures.split(',')]

  h_cols = [t.strip()+" K 1H-RR" for t in temperatures.split(',')]
  n_cols = [t.strip()+" K 15N-RR" for t in temperatures.split(',')]

  # initialize empty dictionaries (fill from results file)
  n_columns = {}
  h_columns = {}

  full_fname = directory + "/" + csvfile
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

    for i in range(len(files)):
      n_columns[h_cols[i].split(" K")[0]] = df[n_cols[i]]
      h_columns[n_cols[i].split(" K")[0]] = df[h_cols[i]]

    ts_int = [float(j) for j in ts]
    for i in range(len(residues)):
      ts_work = []
      n = []
      h = []
      if not(np.isnan(n_columns[n_cols[0].split(" K")[0]][i])) and not(np.isnan(h_columns[h_cols[0].split(" K")[0]][i])):
        for j in range(len(files)):
          if not(np.isnan(n_columns[n_cols[j].split(" K")[0]][i])):
            n.append(n_columns[n_cols[j].split(" K")[0]][i])
          if not(np.isnan(h_columns[h_cols[j].split(" K")[0]][i])):
            h.append(h_columns[h_cols[j].split(" K")[0]][i])
            ts_work.append(ts_int[j])

        # find best linear fit
        ph = np.polyfit(ts_work,h,1)

        # calculate residuals (hres) 
        [f,p, hres] = f_calc(ts_work,h)

        # flag outliers 
        flag = 0
        t_flag = []
        h_flag = []
        for idx, res in enumerate(hres):
          if abs(res) > (2*stddev): # stddev calculated in main_program() 
            t_flag.append(ts_work[idx])
            h_flag.append(h[idx])

        f, (p1,p2) = plt.subplots(2, sharex=True, figsize=(10,10))
        plt.suptitle("Residue "+str(residues[i]))

        # plot 1H chemical shifts vs temperature
        p1.scatter(ts_work,h, s=160, facecolors='none', edgecolors='goldenrod')

        # 1H chemical shifts for review 
        if flag_review == 1:
          p1.scatter(t_flag,h_flag, s=320, facecolors='none', edgecolors='red')

        # plot 1H chemical shift vs. temperature line of best fit
        p1.plot(ts_int,np.polyval(ph, ts_int), color="black", linewidth=2)

        p1.set_ylabel("$\delta$ (ppm)")

        # plot residuals vs temperature
        p2.scatter(ts_work,hres, s=160, facecolors='red', edgecolors='goldenrod')
        p2.set_ylabel("Residual (ppm)")
        p2.set_xlabel("Temperature (K)")

        if numeric_ids:
          chars = len(str(largest_id))
          plt.savefig("residue"+str(residues[i]).zfill(chars)+".png")
        else:
          plt.savefig(str(residues[i]).strip()+".png")

        plt.close()

  return 0

  ## end of plot_figs() function definition


#######################

# f-test to help solve the model selection problem
def f_calc(t,cs):
  lin_fit = np.polyfit(t,cs,1)
  hres = np.polyval(lin_fit, t)-cs
  rss_lin = np.sum((hres)**2)

  quad_fit = np.polyfit(t,cs,2)
  rss_quad = np.sum((np.polyval(quad_fit, t)-cs)**2)

  f = (rss_lin-rss_quad)/(rss_quad/(len(cs)-3))
  p = 1.0 - scipy.stats.f.cdf(f,1,9-2-1)
  return [f,p, hres]

  ## end of f_calc() function definition


########################################
# start of function executed by Celery #
########################################

@app.task
def pyShiftTrack(**kwargs):
  global email
  global jobId
  email = kwargs['email']
  jobId = kwargs['jobId']


  # load config file and set up run
 
  # these definitions are used elsewhere, but don't really need to vary; we hard code them here
  global directory
  global csvfile
  directory= '.'
  csvfile = 'ShiftTrack.csv'
  
  # store current working directory; change to job dir
  start_wd = os.getcwd()
  os.chdir("/home/nmr/shiftt/ShiftTrack/jobs/"+str(jobId))

  # open log
  global log
  log = open("logfile.txt", 'w')

  # load configuration file created during job submission
  global projname
  global assignments
  global peaklists
  global references
  global nom_temperatures
  global start_index
  global n_xi
  global bf_h
  global bf_n
  global dist_cutoff
  global short_line_limit
  global flag_review
  global calc_temp

  config_file = open("STconfig.txt", "rU")
  config = []
  for line in config_file:
    config.append(line.strip())

  if len(config) == 13:
    projname = config[0]
    assignments = config[1]
    peaklists = config[2]
    references = config[3]
    nom_temperatures = config[4]
    start_index = config[5]

    # 15N/1H Ξ (reference compound: liq. NH3)
    n_xi = float(config[6]) # default: 0.10132912
  
    # base transmitter frequencies (Hz)
    bf_h = float(config[7]) # default: 600130000
    bf_n = float(config[8]) # default: 60810645.0

    # distance cut-off, units of ppm squared (see usage below)
    #  - large values will increase runtime and possibly find unlikely 
    #    "solutions" instead of reporting none found
    #  - small values will prevent finding valid solutions with larger
    #    point spacing in the 1H-15N plane
    #  - the sweet spot may depend on both the protein and the ΔT (0.25 has 
    #    been found to work well for Adnectins with ΔT=5 K; 0.125 has been 
    #    found to work well for hisactophilin with ΔT=2.5 K)
    #  - the value supplied via web server submissions is in units of ppm
    #    therefore we square it
    dist_cutoff = float(config[9])**2

    # optionally find lines shorter than the number of temperatures at which we have data...
    # how short is too short?
    short_line_limit = int(config[10])

    flag_review = int(config[11])
    calc_temp = int(config[12])
  else:
    log.write("Aborted due to configuration file error!"+'\n')
    log.close()
    email_error()
    # restore working directory; probably unnecessary
    os.chdir(start_wd)
    return

  # start
  if jobId != "example":
    main_ret = main_program() 
  else:
    main_ret = 0

  # close log 
  log.close()

  # send results
  if main_ret == 0:
    email_results()
  else:
    email_error() 

  # restore working directory; probably unnecessary
  os.chdir(start_wd)
  
  return

