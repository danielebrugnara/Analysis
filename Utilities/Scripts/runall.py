#!/usr/bin/env python3

import os
import re
import io
import sys
import argparse
import glob
import time
import threading

#####Directories of data and configurations
replay_dir=os.getcwd()

data_dir=replay_dir+'/../DataRAW/mugast/e786s'

configs_dir=replay_dir+'/Configs'

#####IP of the grid machines
machines = {
    'agata00':  '10.64.19.19',
    'agata01':  '10.64.19.12',
    'agata02':  '10.64.19.16',
    'agata03':  '10.64.19.17',
    'agata04':  '10.64.19.22',
    'agata05':  '10.64.19.21',
    'agata06':  '10.64.19.9',
    'agata07':  '10.64.19.14',
    'agata08':  '10.64.19.25'
}

Brho = {}

#####Username
username = 'dbrugnara'

#####Matches the run name with the data found in the data directory
def GetDataDirectories(runs):
    directories = {}
    for run in runs:
        for name in glob.glob(data_dir+'/run_*'):
            if 'run_'+run in name:
                directories[run]=name
                print('Found run :', 'run_'+run, ' in directory : ', name)
    if len(directories)!=len(runs):
        print("NOT perfect match between provided runs and data found")
    return directories

#####If asked, removes previous replays folders
def DeletePreviousReplays(runs):
    for run in runs:
        for name in glob.glob(replay_dir+'/Replay_'+run+'*'):
            if 'Replay_'+run in name:
                print("removing : ", name)
                try:
                    os.unlink(name+'/Data')
                except FileNotFoundError:
                    print('No link present in :', name)
                os.system('rm -rf '+name)


#####Sets up new Replay directories
def CreateFolders(data_directories, use_previous):

    replay_directories = {}

    for run in data_directories.keys():
        try:
            directory = replay_dir+'/Replay_'+str(run)
            os.mkdir(directory)
            replay_directories[run]=directory
        except FileExistsError:
            ##In this case the replay is already present, a new directory will be created
            if not use_previous:
                i=1
                while True:
                    directory_added = replay_dir+'/Replay_'+str(run)+'_'+str(i)
                    try:
                        os.mkdir(directory_added)
                        replay_directories[run]=directory_added
                        break
                    except FileExistsError:
                        i=i+1
            else:
                replay_directories[run]=directory
                
    
    ##Copying the configuration files to the replay directory
    configurations=['Conf', 'ConfVAMOS', 'gen_conf.py', 'ADF.conf', 'Topology.conf', 'startReplay.sh', 'Topology_noancillary.conf','startReplay_noancillary.sh']
    if not use_previous:
        for run in replay_directories.keys():
            print('Linking run :', run)
            command = 'ln -s '+ data_directories[run]+'/Data '+replay_directories[run]+'/.'
            os.system(command)
            print('Copying configs of run :', run)
            for conf_file in configurations:
                command = 'cp -r '+configs_dir+'/'+conf_file+' '+replay_directories[run]+'/' 
                os.system(command)
    return replay_directories


#####Runs the genconf on all directories
def RunGenConf(replay_directories):
    for directory in replay_directories.keys():
        print("Running gen_conf for : "+replay_directories[directory])
        command = replay_directories[directory]+'/gen_conf.py'
        os.chdir(replay_directories[directory])
        os.system(command)
        os.chdir(replay_directories[directory])

#####Modify brho based on run number
def ReadBrho():
    ##Read Brho config
    with open(configs_dir+'/BrhoList.txt') as fp:
        line = fp.readline()
        while line:
            values = line.split()
            try:
                if len(values[0])==4:
                    #runs.append(line.strip())
                    brho_value = values[1].strip()
                    if float(brho_value)<2 and float(brho_value)>0.1:
                        print('Setting BRho value for run ', values[0], ' to : ', brho_value)
                        Brho[values[0]] = brho_value
            except IndexError:
                print("Problem parsing Brho config line: ",values)
            line = fp.readline()


def SetBrhos(value, replay_current_dir):
    infile = open(replay_current_dir+'/gen_conf.py')
    outfile = open(replay_current_dir+'/gen_conf_tmp.py', 'w+')
    for line in infile:
        line = re.sub(r'^"BRho\s*[0-9].*', '"Brho\t'+value+'",', line)
        outfile.write(line)
    infile.close()
    outfile.close()
    os.rename(replay_current_dir+'/gen_conf.py',replay_current_dir+'/gen_conf_old.py' )
    os.rename(replay_current_dir+'/gen_conf_tmp.py',replay_current_dir+'/gen_conf.py' )
    os.system('chmod +x '+replay_current_dir+'/gen_conf.py')

#####Sets up threads and allocates work on the machines
def main():
    #####Setting up the parser
    parser = argparse.ArgumentParser(description='Script to run replays in Runs.txt')
    parser.add_argument('--cpus',
                    type=int,
                    nargs='+',
                    dest='cpus',
                    required=True,
                    help='agata machines to run replay on')
    parser.add_argument('--delete_previous',
                    action='store_true',
                    dest='delete_previous',
                    help='deletes all replays of a given run number')
    parser.add_argument('--use_previous',
                    action='store_true',
                    dest='use_previous',
                    help='use previous replay of type Replay_****, keeps the old configs')

    args = parser.parse_args()

    #####File containing the runs to run the replay on
    filepath = configs_dir+"/Runs.txt"
    #####Reading the requested runs
    runs = []
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            if len(line)==5:
                runs.append(line.strip())
            line = fp.readline()

    for run in runs:
        print('Loaded runs :', run)

    #####Mathcing requests with data folders
    data_directories = GetDataDirectories(runs)

    #####Delete previous replays?
    if (args.delete_previous and not args.use_previous):
        DeletePreviousReplays(runs)

    #####Creating folders
    replay_directories = CreateFolders(data_directories, args.use_previous)
    
    
    #####Running gen conf and reading Brho
    if not args.use_previous:
        ReadBrho()
        for run in runs:
            try:
                SetBrhos( Brho[run], replay_directories[run])
            except:
                print("Unable to find run in dictionaries (no brho in dictionary for this run", run,"), running with default: ")
        RunGenConf(replay_directories)


    #####Stack of jobs to complete
    stack = []
    for run in runs:
        stack.append(run)

    #####List of available threads
    threads = list()
    #####Lock instantiation to access common resource (stack of jobs)
    lock = threading.Lock()

    #####Thread friendly access to shared resource
    def GetJob():
        lock.acquire()
        try:
            run = stack.pop()
        finally:
            ##Releasing the lock once the stack is popped
            lock.release()
        return run

    #####If some crysyals are not on the data, a new topology will be created, excluding them
    def RemoveFromTopology(crystal, replay_current_dir):
        infile = open(replay_current_dir+'/Topology.conf')
        outfile = open(replay_current_dir+'/Topology_tmp.conf', 'w+')
        for line in infile:
            line = line.replace(str(crystal+" "), "")
            outfile.write(line)
        infile.close()
        outfile.close()
        os.rename(replay_current_dir+'/Topology.conf',replay_current_dir+'/Topology_old.conf' )
        os.rename(replay_current_dir+'/Topology_tmp.conf',replay_current_dir+'/Topology.conf' )

    def TouchPSAFile(crystal, replay_current_dir):
        command = "touch "+replay_current_dir+"/Data/"+crystal+"/SRM_AGATA_psa_0000.adf"
        os.system(command) 
 
    #####Is the current job in the runs to exclude in the stop file?
    def Stop(current_run):
        with open(configs_dir+'/StopFile.txt', 'r') as stop_file:
            line = stop_file.readline()
            while line:
                print(line)
                if len(line)==5:
                    if line.strip()==current_run:
                        print("Stopping run :", current_run)
                        return True
                    else:
                        line = stop_file.readline()
            return False

    ######SSH call for a job
    def Execute(machine_nr):
        ##Access to shared resource
        current_run = GetJob()
        
        replay_current_dir = replay_directories[current_run]
        command = 'ssh '+username+'@'+machines['agata0'+str(machine_nr)]+' \'cd '+replay_current_dir+'; sh startReplay.sh\''
        command_noancillary = 'ssh '+username+'@'+machines['agata0'+str(machine_nr)]+' \'cd '+replay_current_dir+'; sh startReplay_noancillary.sh\''

        ##Removing previous logs
        for name in glob.glob(replay_current_dir+'/*log.txt'):
            print ("Removing logs:", name)
            try:
                os.remove(name)
            except OSError as error:
    	        print(error) 
    	        print("File path can not be removed") 
        
        
        ##Loop to identify missing crstals
        iteration = 0 
        while True:
            logs = open(replay_current_dir+"/nr_"+str(iteration)+"_log.txt", "w+")
            output = os.popen(command).read()
            logs.write(output)
            logs.close()
            tailcommand = 'tail -20 '+logs.name
            print("printing output from :", tailcommand)
            lines = os.popen(tailcommand).readlines()
            if Stop(current_run):
                break

            ##Error handling of the launched command through cout redirection
            if '::process_initialise' in lines[18]:
                ##Error, missing data for a given crystal (likely)
                error_log = open(replay_current_dir+'/ERROR_log.txt', "a+")
                error = lines[18].split('/Data/')
                error_log.write('------->Unable to find crystal : ')
                missing_crystal = error[1][:3]
                error_log.write(missing_crystal)
                error_log.write('\n')
                #error_log.write(lines[18])
                error_log.close()
                if missing_crystal == "vam":
                    command = command_noancillary
                    os.system("touch "+replay_current_dir+"/NOANCILLARYINRUN.txt")
                else:
                    print("Touching missing psa.adf for crystal:", missing_crystal)
                    TouchPSAFile(missing_crystal, replay_current_dir)
                ##Removing missing crystal from topology and trying again
                
                #the following will remove it completely from the topology
                #print("Removing missing crystal from topology :", missing_crystal)
		#RemoveFromTopology(missing_crystal, replay_current_dir)
		
                iteration+=1
                if iteration>30:
                    break
                continue
            if 'TreeBuilder' in lines[19]:
                ##All is ok
                break
            else:
                ##Unknown exception, writing to ERROR_log.txt
                error_log = open(replay_current_dir+'/ERROR_log.txt', "a+")
                error_log.write("------->Unrecognized error:\n")
                for line in lines:
                    error_log.write(line)
                error_log.close()
                break

    #####Work untill stack is empty
    def ThreadWorker(machine_nr):
        try:
            while True:
                Execute(machine_nr)
        except IndexError:
            print("Finished runs")
        return
        
    #####Cpus available from parser
    for cpu in args.cpus:
        if 'agata0'+str(cpu) in machines.keys():
            x = threading.Thread(target=ThreadWorker, args=(cpu,))
            threads.append(x)
        else:
            print("Provided machine number not existent, skipping cpu number : ", cpu)
        
    #####Starting all threads
    for thread in threads:
        thread.start()

    #####Waiting for alla threads to finish
    for thread in threads:
        thread.join()
    
    print('All threads have joined correctly, you may now go in peace')

if __name__ == '__main__':
    main()
