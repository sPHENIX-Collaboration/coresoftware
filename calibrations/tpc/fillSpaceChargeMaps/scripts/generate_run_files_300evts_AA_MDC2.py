#!/usr/bin/env python

# coding=utf-8

from pathlib import Path
import stat

def chmod_px(name):
    f = Path(name)
    f.chmod(f.stat().st_mode | stat.S_IEXEC)
    return 0

work_dir = "/sphenix/user/shulga/Work/TpcPadPlane_phi_coresoftware/coresoftware/calibrations/tpc/fillSpaceChargeMaps/"
init_str = "Initialdir      = "+work_dir

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
    init_str,
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

ff= open("./run_all_300evts_AA_jobs_MDC2.sh","w+")
ff.write("#!/usr/bin/bash"+"\n")
ff_all_dag= open("./run_all_300evts_AA_jobs_MDC2_dag.sh","w+")
ff_all_dag.write("#!/usr/bin/bash"+"\n")
evt_start = [0,  75, 155, 233, 312, 392, 471, 551, 630, 712, 793, 873, 953, 1033, 1113, 1194, 1272, 1352, 1432, 1511, 1592, 1670, 1750, 1830, 1910, 1990, 2069, 2149, 2230, 2310, 2391, 2470, 2550, 2631, 2711, 2792, 2871, 2952, 3032, 3113, 3193, 3273, 3352, 3433, 3514, 3595, 3676, 3755, 3834, 3915, 3995, 4076, 4156, 4235, 4315, 4396, 4476, 4557, 4636, 4718, 4798, 4877, 4957, 5037, 5118, 5200, 5278, 5359, 5439, 5518, 5598, 5680, 5760, 5841, 5922, 6002, 6081, 6162, 6241, 6321, 6400, 6481, 6561, 6641, 6722, 6804, 6884, 6964, 7045, 7126, 7206, 7286, 7366, 7447, 7527, 7606, 7687, 7767, 7847, 7928, 8008, 8087, 8167, 8248, 8328, 8408, 8489, 8566, 8646, 8727, 8809, 8890, 8970, 9049, 9129, 9209, 9288, 9368, 9449, 9529, 9608, 9688, 9769, 9849]
evt_end = [ 170, 248, 327, 407, 486, 566, 645, 727, 808, 888, 968, 1048, 1128, 1209, 1287, 1367, 1447, 1526, 1607, 1685, 1765, 1845, 1925, 2005, 2084, 2164, 2245, 2325, 2406, 2485, 2565, 2646, 2726, 2807, 2886, 2967, 3047, 3128, 3208, 3288, 3367, 3448, 3529, 3610, 3691, 3770, 3849, 3930, 4010, 4091, 4171, 4250, 4330, 4411, 4491, 4572, 4651, 4733, 4813, 4892, 4972, 5052, 5133, 5215, 5293, 5374, 5454, 5533, 5613, 5695, 5775, 5856, 5937, 6017, 6096, 6177, 6256, 6336, 6415, 6496, 6576, 6656, 6737, 6819, 6899, 6979, 7060, 7141, 7221, 7301, 7381, 7462, 7542, 7621, 7702, 7782, 7862, 7943, 8023, 8102, 8182, 8263, 8343, 8423, 8504, 8581, 8661, 8742, 8824, 8905, 8985, 9064, 9144, 9224, 9303, 9383, 9464, 9544, 9623, 9703, 9784, 9864, 9943, 10000]
evt_bX =  [1508071.0, 3016509.0, 4524020.0, 6032112.0, 7540028.0, 9048092.0, 10556072.0, 12064371.0, 13572143.0, 15080178.0, 16588072.0, 18096105.0, 19604092.0, 21112166.0, 22620146.0, 24128151.0, 25636093.0, 27144133.0, 28652094.0, 30160125.0, 31668120.0, 33176312.0, 34684455.0, 36192201.0, 37700299.0, 39208338.0, 40716228.0, 42224205.0, 43732346.0, 45240219.0, 46748391.0, 48256211.0, 49764321.0, 51272217.0, 52780276.0, 54288364.0, 55796155.0, 57304180.0, 58812172.0, 60320612.0, 61828166.0, 63336720.0, 64844207.0, 66352327.0, 67860452.0, 69368220.0, 70876499.0, 72384379.0, 73892278.0, 75400246.0, 76908483.0, 78416343.0, 79924240.0, 81432250.0, 82940306.0, 84449218.0, 85956368.0, 87464432.0, 88972282.0, 90480292.0, 91988279.0, 93496320.0, 95004286.0, 96512289.0, 98020429.0, 99528306.0, 101036272.0, 102544367.0, 104052284.0, 105560289.0, 107068427.0, 108576408.0, 110084415.0, 111592686.0, 113100369.0, 114608368.0, 116116420.0, 117624429.0, 119132375.0, 120640371.0, 122148359.0, 123656423.0, 125164763.0, 126673264.0, 128180360.0, 129688558.0, 131196596.0, 132704683.0, 134212939.0, 135720446.0, 137228683.0, 138736381.0, 140244545.0, 141752581.0, 143260466.0, 144768499.0, 146276523.0, 147784726.0, 149292502.0, 150800480.0, 152308621.0, 153816539.0, 155324424.0, 156832627.0, 158340668.0, 159848664.0, 161356504.0, 162864630.0, 164372489.0, 165880652.0, 167388456.0, 168896455.0, 170404493.0, 171912555.0, 173421358.0, 174928568.0, 176436571.0, 177945044.0, 179452976.0, 180961051.0, 182468742.0, 183976536.0, 185484542.0, 186992509.0]

for j, (start,end) in enumerate(zip(evt_start,evt_end)):
    filename_bx = "./condor_macros/run_all_300evts_AA_jobs_MDC2_{}.sh".format(int(evt_bX[j]))
    filename_dag = "./condor_macros/run_all_300evts_AA_jobs_MDC2_{}.dag".format(int(evt_bX[j]))
    ff_bx= open(filename_bx,"w+")
    ff_dag= open(filename_dag,"w+")
    ff_bx.write("#!/usr/bin/bash"+"\n")
    ff.write("source {}\n".format(filename_bx))
    ff_all_dag.write("condor_submit_dag  {}\n".format(filename_dag))
    parent_str = "PARENT"
    for i in range(start,end+10):
        filename = "./condor_macros/run_files_300evts_AA_MDC2_{}_{}.sh".format(j,i)
        f= open(filename,"w+")
        f.write("#!/usr/bin/bash"+"\n")
        f.write("source macros/run_files_300evts_AA_MDC2.sh {} {} {}".format(i,i+1,evt_bX[j])+"\n")
        f.close
        chmod_px(filename)
        filename_job = "./condor_macros/condor_run_files_300evts_AA_MDC2_{}_{}.job".format(j,i)
        ff_bx.write("condor_submit {}".format(filename_job)+"\n")
        ff_dag.write("JOB A_{} {}\n".format(i,filename_job ))
        parent_str += " A_{}".format(i)
        f_job= open(filename_job,"w+")
        n_line = 0
        for lines in introduction:
            f_job.write(lines+"\n")
            if n_line==3:
                f_job.write("# The executable we want to run."+"\n")
                f_job.write("Executable      = condor_macros/run_files_300evts_AA_MDC2_{}_{}.sh".format(j,i)+"\n")
                f_job.write(""+"\n")
                f_job.write(""+"\n")
                f_job.write("# The argument to pass to the executable."+"\n")
                f_job.write("Arguments       = \"run job 300 evts AA {} {}\"".format(j,i)+"\n")
            if n_line==38:
                f_job.write("# The job's stdout is sent to this file."+"\n")
                f_job.write("Output          = " + work_dir + "Out/myjob_300evts_AA_{}_{}.out".format(j,i)+"\n")
                f_job.write(""+"\n")
                f_job.write("# The job's stderr is sent to this file."+"\n")
                f_job.write("Error           = " + work_dir + "Out/myjob_300evts_AA_{}_{}.err".format(j,i)+"\n")
                f_job.write(""+"\n")
                f_job.write("# The condor log file for this job, useful when debugging."+"\n")
                f_job.write("Log             = " + work_dir + "Out/condor_300evts_AA_{}_{}.log".format(j,i)+"\n")
                f_job.write(""+"\n")
    
            n_line+=1

        f_job.close
    ff_bx.close
    chmod_px(filename_bx)
    
    #chmod_px(filename)
    filename_B = "./condor_macros/add_hist_{}.job".format(int(evt_bX[j]))
    ff_dag.write("JOB B {}\n".format(filename_B ))
    ff_dag.write(parent_str + " CHILD B")
    ff_dag.close
    chmod_px(filename_dag)

    # files for adding histos
    # - .sh:
    filename_B_sh = "./condor_macros/add_hist_{}.sh".format(int(evt_bX[j]))
    ff_B_sh = open(filename_B_sh,"w+")

    ff_B_sh.write("#!/usr/bin/bash"+"\n")
    ff_B_sh.write("source /opt/sphenix/core/bin/sphenix_setup.sh -n new"+"\n")
    ff_B_sh.write("./macros/add_histos_bX.py {} {} \n".format(j,int(evt_bX[j])))

    ff_B_sh.close
    chmod_px(filename_B_sh)

    # - .job:
    add_file_content = [
        "Universe        = vanilla",
        "Executable      = {}".format(filename_B_sh),
        "Arguments       = \"run job add histos AA {}\"".format(int(evt_bX[j])),
        "request_memory = 7.1GB",
        init_str,
        "# The job's stdout is sent to this file.",
        "Output          = ./Out/myjob_add_histos_AA_{}.out".format(int(evt_bX[j])),
        "# The job's stderr is sent to this file.",
        "Error           = ./Out/myjob_add_histos_AA_{}.err".format(int(evt_bX[j])),
        "# The condor log file for this job, useful when debugging.",
        "Log             = ./Out/condor_add_histos_AA_{}.log".format(int(evt_bX[j])),
        "Queue"
        ]
    ff_B = open(filename_B,"w+")
    for line in add_file_content:
        ff_B.write(line+"\n")
    ff_B.close


ff.close
ff_all_dag.close
