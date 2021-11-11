from subprocess import Popen, PIPE, STDOUT
import os.path
import pandas as pd
import io
#import MDAnalysis as mda
#from MDAnalysis.analysis import align

#  constants for gromacs commands
GMX_CALPHA='3'
GMX_PROTEIN='1'
GMX_BACKBONE='4'


## UTILITIES

# generic commands
def run_cmd(cmd,options):
    p = Popen(cmd, stdin=PIPE, stderr=PIPE)
    if (options != ''):
        stderr = p.communicate(input=options.encode())[0]

# test if files exists
def test_files(filelist):
    for filename in filelist:
        if not (os.path.isfile(filename)):
            print("File "+filename+" not found")

## GROMACS specific routines
def read_xvg(xvg,comment='#',comment2="@",sep='\s+'):
    ''' Reads an xvg file created by Gromacs and outputs of a Pandas dataframe 
        skipping the comment lines. The calling process should check the column names
        (default x and y)
    '''
    lines = "".join([line for line in open(xvg) 
                     if not (line.startswith(comment) or line.startswith(comment2))])
    return pd.read_csv(io.StringIO(lines), sep=sep, header=0,names=['x','y'])


def gmx_covar(xtc,tpr,start,eigenvec,options):

    cmd=['gmx','covar','-f',xtc,'-s',tpr,'-b', start,'-v',eigenvec]
    run_cmd(cmd,options)
    

def gmx_anaeig(xtc,tpr,start,eigenvec,options,xvg='FirstPlane.xvg'):

    cmd=['gmx','anaeig','-f',xtc,'-s',tpr ,'-b', start, '-2d',xvg,'-first','1','-last','2','-v',eigenvec]
    run_cmd(cmd,options)
    df=pd.read_csv(xvg,delimiter='\s+', skiprows=17,names=['x','y'])
    return df

def gmx_cluster(xtc,tpr,index,start,skip,clindex,options):
    cutoff='0.2'
    cmd=['gmx','cluster','-f',xtc,'-s',tpr ,'-b', start,'-g','-dist','-ev','-sz','-cl','-method', \
     'gromos','-cutoff',cutoff,'-wcl','3','-clndx',clindex,'-skip',skip] # ['-n',index] for holo
    run_cmd(cmd,options)
    df=pd.read_csv('clust-size.xvg',delimiter='\s+', skiprows=18,names=['x','y'])
    return df

def gmx_extract_cluster(xtc,tpr,index,clindex,prefix,options):
    cmd=['gmx','extract-cluster','-f',xtc,'-s',tpr ,'-n',index,'-clusters',clindex,'-o',prefix]
    run_cmd(cmd,options)

def gmx_rmsf(xtc,tpr,rmsf,options):
    cmd=['gmx','rmsf','-f',xtc ,'-s',tpr ,'-res', '-o',rmsf]
    run_cmd(cmd,options)
    df=pd.read_csv(rmsf,delimiter='\s+', skiprows=18,names=['x','y'])
    return df

def get_cluster(cluster_file='cluster.log',max_clusters=5):
    ''' Reads a GROMACS cluster.log file and returns a list of times, one for each cluster'''

    df=pd.read_csv(cluster_file,skiprows=12,delimiter='|')
    df.rename(columns=lambda x: x.strip(),inplace=True)

    # remove spaces from cluster name column
    df['cl.']=df['cl.'].str.strip()

    # set the cluster name for all rows
    for i in range(len(df)):
    	clstr=df.loc[i,'cl.']
    	if (clstr !=''):
            cl=clstr
    	else:
            df.loc[i,'cl.']=cl

    dfg=df.groupby('cl.') # group by cluster names
    cluster_times=[]

    for cluster in range(1,max_clusters+1):
        g=dfg.get_group(str(cluster))['cluster members']
        pts=[]
        # the following splits up each line and adds the times to the list for each cluster
        [[pts.append(float(pt)) for pt in line.split()] for line in g]
        cluster_times.append(pts)
    return (cluster_times)



## MDAnalysis routines

def rmsf(u,pdb):
    '''Calculates RMSF using MDanalysis. Also possible in parallel with PMDA'''
    protein = u.select_atoms("protein")

    # TODO: Need to center and make whole (this test trajectory
    # contains the protein being split across periodic boundaries
    # and the results will be WRONG!)

    # Fit to the initial frame to get a better average structure
    # (the trajectory is changed in memory)
    prealigner = align.AlignTraj(u, u, select="protein and name CA",
                                         in_memory=True).run()
    # ref = average structure
    ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
    # Make a reference structure (need to reshape into a
    # 1-frame "trajectory").
    ref = mda.Merge(protein).load_new(ref_coordinates[:, None, :],order="afc")
    aligner = align.AlignTraj(u, ref, select="protein and name CA",
                                      in_memory=True).run()
    # need to write the trajectory to disk for PMDA 0.3.0 (see issue #15)
    with mda.Writer("rmsfit.xtc", n_atoms=u.atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(u.atoms)
    u = mda.Universe(pdb, "rmsfit.xtc")
    calphas = protein.select_atoms("protein and name CA")

    #rmsfer = rms.RMSF(calphas).run()
    return calphas

