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
evt_start = [0,  75, 155, 233, 312, 392, 471, 551, 630, 712, 793, 873, 953, 1033, 1113, 1194, 1272, 1352, 1432, 1511, 1592, 1670, 1750, 1830, 1910, 1990, 2069, 2149, 2230, 2310, 2391, 2470, 2550, 2631, 2711, 2792, 2871, 2952, 3032, 3113, 3193, 3273, 3352, 3433, 3514, 3595, 3676, 3755, 3834, 3915, 3995, 4076, 4156, 4235, 4315, 4396, 4476, 4557, 4636, 4718, 4798, 4877, 4957, 5037, 5118, 5200, 5278, 5359, 5439, 5518, 5598, 5680, 5760, 5841, 5922, 6002, 6081, 6162, 6241, 6321, 6400, 6481, 6561, 6641, 6722, 6804, 6884, 6964, 7045, 7126, 7206, 7286, 7366, 7447, 7527, 7606, 7687, 7767, 7847, 7928, 8008, 8087, 8167, 8248, 8328, 8408, 8489, 8566, 8646, 8727, 8809, 8890, 8970, 9049, 9129, 9209, 9288, 9368, 9449, 9529, 9608, 9688, 9769, 9849]
evt_end = [ 170, 248, 327, 407, 486, 566, 645, 727, 808, 888, 968, 1048, 1128, 1209, 1287, 1367, 1447, 1526, 1607, 1685, 1765, 1845, 1925, 2005, 2084, 2164, 2245, 2325, 2406, 2485, 2565, 2646, 2726, 2807, 2886, 2967, 3047, 3128, 3208, 3288, 3367, 3448, 3529, 3610, 3691, 3770, 3849, 3930, 4010, 4091, 4171, 4250, 4330, 4411, 4491, 4572, 4651, 4733, 4813, 4892, 4972, 5052, 5133, 5215, 5293, 5374, 5454, 5533, 5613, 5695, 5775, 5856, 5937, 6017, 6096, 6177, 6256, 6336, 6415, 6496, 6576, 6656, 6737, 6819, 6899, 6979, 7060, 7141, 7221, 7301, 7381, 7462, 7542, 7621, 7702, 7782, 7862, 7943, 8023, 8102, 8182, 8263, 8343, 8423, 8504, 8581, 8661, 8742, 8824, 8905, 8985, 9064, 9144, 9224, 9303, 9383, 9464, 9544, 9623, 9703, 9784, 9864, 9943, 10000]
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
