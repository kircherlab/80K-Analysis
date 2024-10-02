"""
Description:            This script extracts PWMs from e.g. JASPAR data base to pass them to MEME FIMO in a subsequent step

Output:                 Motif PWMs of interest from JASPAR

example bash command:   python ./allMotifsWithSignificantEffects/removeSimilarMotifs.py --outputPWM condensed_test.txt
"""


import click
# import h5py
import numpy as np
import os
import pandas as pd

@click.command()
@click.option(
    "--PWM_file",
    "PWM_file",
    required=True,
    multiple=False,
    type=str,
    default="test.txt",
    help="e.g. meme file",
)
@click.option(
    "--outputPWM",
    "out_pwm",
    required=True,
    multiple=False,
    type=str,
    help="output PWM file",
)
@click.option(
    "--tomtomSig",
    "tomtom_qual_level",
    required=True,
    multiple=False,
    default=0.01,
    help="alpha for sig testing for tomtom matches",
)
def cli(PWM_file, out_pwm, tomtom_qual_level):
    # read PWM file
    data_base=open(PWM_file, "r")
    data_base_entries=data_base.read().split("MOTIF")
    data_base_entries=data_base_entries[1:] #erster entry leer
    data_base.close()

    #initialize outfile with significant PWMs in outPWM  set up
    outPWM=open(out_pwm, "w")
    background = 0.25

    outPWM.write('MEME version 4\n\n')
    outPWM.write('ALPHABET= ACGT\n\n')
    outPWM.write('strands: + -\n\n')
    outPWM.write('Background letter frequencies (from unknown source):\n')
    outPWM.write('A %.3f C %.3f G %.3f T %.3f\n\n' % (background, background, background, background))

    # look for particular motif and write them to output file
    var_len_data_base_entries=len(data_base_entries)
    p=0
    while p < var_len_data_base_entries:
        #write entry p to meme
        # create temp meme file for PWM to be tested using tomtom againgst rest of PWMS
        tmp_file=open("tmp.txt", "w")
        background = 0.25
        tmp_file.write('MEME version 4\n\n')
        tmp_file.write('ALPHABET= ACGT\n\n')
        tmp_file.write('strands: + -\n\n')
        tmp_file.write('Background letter frequencies (from unknown source):\n')
        tmp_file.write('A %.3f C %.3f G %.3f T %.3f\n\n' % (background, background, background, background))
        append_to_meme_file(data_base_entries[p], tmp_file)
        tmp_file.close()



        # run tomtom of motif p aiannst rest
        cmd = 'export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.4:$PATH'
        os.system(cmd)
        cmd = 'tomtom -no-ssc -oc . --verbosity 1 -text -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 tmp.txt %s > %s' % (PWM_file , "tmp_tomtom")
        os.system(cmd)

        # select good tomtom matches with qvalue < tomtom_qual_level
        if os.path.getsize("tmp_tomtom")!=0:
            df_tomtom=pd.read_csv("tmp_tomtom", sep="\t", skipfooter=4)
            #print(df_tomtom["q-value"].shape)
            df_tomtom_sig=df_tomtom.loc[df_tomtom['q-value'] < tomtom_qual_level]
            #print(df_tomtom_sig["q-value"].shape)

            # remove similar motifs
            similarMotifs=[]
            for index, row in df_tomtom_sig.iterrows():
                tmp_motif=row['Target_ID']
                #print("target", tmp_motif)
                #print("querey", )
                if row['Target_ID']!=row['Query_ID']:
                    for g in range (p, len(data_base_entries), 1):
                        #print(data_base_entries[g].split("\n")[0].split()[0])
                        if data_base_entries[g].split("\n")[0].split()[0] == str(tmp_motif):
                            similarMotifs.append(data_base_entries[g])
                            #print(data_base_entries[g])


            #data_base_entries = list(set(data_base_entries).difference(set(similarMotifs)))[:]
            c = [x for x in data_base_entries if not x in similarMotifs]
            data_base_entries=c[:]
            #print (len(similarMotifs), "motifs removed")
        #print ("list of removed motifs", similarMotifs)
        var_len_data_base_entries=len(data_base_entries)
        p=p+1

    # write remaining PWMS to out meme file
    print("remaining motifs", len(data_base_entries))
    for g in range (0, len(data_base_entries), 1):
        #print(data_base_entries[g])
        append_to_meme_file(data_base_entries[g], outPWM)
    outPWM.close()



def append_to_meme_file(data_base_entry, existingFileName):
    #print("MOTIF " + data_base_entry)
    existingFileName.write("MOTIF " + data_base_entry)

if __name__ == "__main__":
    cli()
