#!/usr/bin/env python

# coding=utf-8

introduction = [
    "# All local jobs are part of the vanilla universe.",
    "Universe        = vanilla",
    "",
    "# The requirement line specifies which machines we want to",
    "# run this job on.  Any arbitrary classad expression can",
    "# be used.",
    "#Requirements    = (CPU_Speed >= 1)",
    "",
    "# Rank is an expression that states how to rank machines which ",
    "# have already met the requirements expression.  Essentially, ",
    "# rank expresses preference.  A higher numeric value equals better ",
    "# rank.  Condor will give the job the machine with the highest rank.",
    "#    Rank = CPU_Speed",
    "",
    "# Jobs by default get 1.4Gb of RAM allocated, ask for more if needed",
    "# but if a job needs more than 2Gb it will not be able to run on the",
    "# older nodes",
    "request_memory = 7.1GB",
    "",
    "# If you need multiple cores you can ask for them, but the scheduling",
    "# may take longer the \"larger\" a job you ask for",
    "request_cpus = 1",
    "",
    "# This flag is used to order only one's own submitted jobs ",
    "# The jobs with the highest numbers get considered for ",
    "# scheduling first.",
    "#Priority        = 4",
    "",
    "# Copy all of the user's current shell environment variables ",
    "# at the time of job submission.",
    "#GetEnv          = True",
    "",
    "# Used to give jobs a directory with respect to file input ",
    "# and output.",
    "Initialdir      = /sphenix/user/shulga/Work/IBF/readDigitalCurrents/",
    "",
    "# Input file given to the job.",
    "#Input           = /dev/null",
    "",
    "",
    "# This should be the last command and tells condor to queue the",
    "# job.  If a number is placed after the command (i.e. Queue 15)",
    "# then the job will be submitted N times.  Use the $(Process)",
    "# macro to make your input/output and log files unique.",
    "Queue"
]

ff= open("./run_all_AA_jobs.sh","w+")
ff.write("#!/usr/bin/bash"+"\n"),
#evt_start = [0,8,16,23,31,39,47,55,63,71,79,87]
#evt_end = [9,17,24,32,40,48,56,64,72,80,88,96]
evt_start = [0,80,160,230,310,390,470,550,630,710,790,870]
evt_end = [90,170,240,320,400,480,560,640,720,800,880,960]
evt_bX = [1508071.0, 3016509.0, 4524020.0, 6032112.0, 7540028.0, 9048092.0, 10556072.0, 12064371.0, 13572143.0, 15080178.0, 16588072.0, 18096105.0]
for j, (start,end) in enumerate(zip(evt_start,evt_end)):
    for i in range(start,end):
        filename = "./run_macros/run_files_AA_{}_{}.sh".format(j,i)
        f= open(filename,"w+")
        f.write("#!/usr/bin/bash"+"\n")
        f.write("source macros/run_files_AA.sh {} {} {}".format(i,i+1,evt_bX[j])+"\n")
        f.close
        filename_job = "./run_macros/condor_run_files_AA_{}_{}.job".format(j,i)
        ff.write("condor_submit {}".format(filename_job)+"\n")
        f_job= open(filename_job,"w+")
        n_line = 0
        for lines in introduction:
            f_job.write(lines+"\n")
            if n_line==3:
                f_job.write("# The executable we want to run."+"\n")
                f_job.write("Executable      = run_macros/run_files_AA_{}_{}.sh".format(j,i)+"\n")
                f_job.write(""+"\n")
                f_job.write(""+"\n")
                f_job.write("# The argument to pass to the executable."+"\n")
                f_job.write("Arguments       = \"run DST HISTO job AA {} {}\"".format(j,i)+"\n")
            if n_line==38:
                f_job.write("# The job's stdout is sent to this file."+"\n")
                f_job.write("Output          = /sphenix/user/shulga/Work/IBF/readDigitalCurrents/Out/myjob_AA_{}_{}.out".format(j,i)+"\n")
                f_job.write(""+"\n")
                f_job.write("# The job's stderr is sent to this file."+"\n")
                f_job.write("Error           = /sphenix/user/shulga/Work/IBF/readDigitalCurrents/Out/myjob_AA_{}_{}.err".format(j,i)+"\n")
                f_job.write(""+"\n")
                f_job.write("# The condor log file for this job, useful when debugging."+"\n")
                f_job.write("Log             = /sphenix/user/shulga/Work/IBF/readDigitalCurrents/Out/condor_AA_{}_{}.log".format(j,i)+"\n")
                f_job.write(""+"\n")
    
            n_line+=1
        f_job.close
ff.close
