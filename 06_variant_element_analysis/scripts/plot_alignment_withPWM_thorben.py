"""
Description:


Example commands:   python VariantEffects/plot_alignment_withPWM.py
Outputs:            fasta sequence file with sequences cooresponding to variants with top x % of all variant effects. CSV with haplo seqs and variant effect

"""

#set seed and load dependencies
import click
import logomaker
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

col_alt_id = 'tmp_variant_fimo_hit_ALT_ID'
col_ref_id = 'tmp_variant_fimo_hit_REF_ID'
col_variant_pos = 'tmp_variant_fimo_hit_variant_pos'
col_motif_id = 'tmp_fimo_results_motif_id'
col_start = 'tmp_fimo_results_start'
col_stop = 'tmp_fimo_results_stop'
col_strand = 'tmp_fimo_results_strand'
col_fimo_q_value = 'tmp_fimo_results_q-value'
col_fimo_matched_sequence = 'tmp_fimo_results_matched_sequence'
col_ref_sequence = 'ref_sequence'
col_alt_sequence = 'alt_sequence'
col_mpralm_logFC = 'tmp_variant_BECALM_logFC'

# TODO: Add strand to ID:
# TODO: move plot to the right to fit the motif better
# TODO: Make motifs clickable to hocomoco (directly putting a html link into the svg using the y_label is not possible)
@click.command()
@click.option(
    "--aligMotifFile",
    "alig_motif_file",
    required=True,
    multiple=False,
    type=str,
    default="results/diffTeloHEAC_CTRL_vs_6h/MPRAlm_Alec_diffTeloHEAC_CTRL_vs_6h_withMotifs.csv",
    help="e.g. starrseq-all-final-toorder_oligocomposition.csv (output of FIMO)",
)
@click.option(
    "--PWMsFile",
    "PWM_file",
    required=True,
    multiple=False,
    type=str,
    default="results/diffTeloHEAC_CTRL_vs_6h/TeloHAEC_CTRL_sigPWMs.txt",
    help="e.g. 2023-01-10_22-29-33 myCounts.minDNAfilt.depthNorm.keepHaps - starr.haplotypes.oligo1.txt (the motif file)",
)
@click.option(
    "--chunk_size",
    "chunk_size",
    required=True,
    multiple=False,
    type=int,
    default=5,
    help="Create subfigures for faster plotting (default: 5)",
)
@click.option(
    "--plottingBackground",
    "plotting_background",
    required=False,
    multiple=False,
    type=float,
    default=0.25,
    help="The plotting background (default: 0.25) will not plot any letters with a frequency below the given value in the sequence logo",
)
@click.option(
    "--variantEffectColumn",
    "col_mpralm_logFC",
    required=False,
    multiple=False,
    type=str,
    default="",
    help="If set the given table is expected to have a column with the given name and the logFC values from e.g. MPRAlm are given for each variant",
)
@click.option(
    "--hocomoco",
    "hocomoco",
    required=False,
    multiple=False,
    type=bool,
    default=False,
    help="If hocomoco is set to true the motifs will be linked to hocomoco website (not working yet)",
)
@click.option(
    "--positiveVariantColor",
    "positive_variant_color",
    required=False,
    type=str,
    default="cyan",
    help="If VariantEffectColumn is set the color of the positive variant effect if no variantEffectColumn is given it is the default color of the variant background",
)
@click.option(
    "--negativeVariantColor",
    "negative_variant_color",
    required=False,
    type=str,
    default="m",
    help="If VariantEffectColumn is set the color of the negative variant effect",
)
@click.option(
    "--variantPosZeroBased",
    "variant_pos_zero_based",
    required=False,
    default=False,
    type=bool,
    help="Use if the given variant position is 0 based (recommended) (default: 1 based because of backward compatibility)",
)
@click.option(
    "--output",
    "out",
    required=True,
    multiple=False,
    type=click.Path(writable=True),
    default="This is a required field",
    help="output name (<chunk number>.svg is added)",
)

def cli(alig_motif_file, PWM_file, chunk_size, col_mpralm_logFC, hocomoco, positive_variant_color, negative_variant_color, plotting_background, variant_pos_zero_based, out):

    #import labels; abs diff activities of seq with ID1 - ID2
    df_aligned=pd.read_csv(alig_motif_file, sep="\t", low_memory=False)
    print("number of variant pairs :", df_aligned.shape[0])
    _text_size=15

    data_base=open(PWM_file, "r")
    data_base_entries=data_base.read().split("MOTIF")
    data_base.close()

    # batch variant motif matches
    chunkNumber=1
    for chunk in np.array_split(df_aligned, chunk_size):
        fig, axs= False, False
        chunk=chunk.reset_index(drop=True)
        chunkNumber=chunkNumber+1
        height_ratio_list=[]
        for i in range (0, chunk.shape[0]*3, 3):
            height_ratio_list.append(3)
            height_ratio_list.append(0.5)
            height_ratio_list.append(0.5)
        fig, axs = plt.subplots(nrows= chunk.shape[0]*3, ncols=1, sharex=True, sharey=False, height_ratios=height_ratio_list, figsize=(16, 3*chunk.shape[0] ))
        if chunk.shape[0]==1: # bug repair, da sonst bei nur einem motif ax[i] mit i=0 nicht gefunden wird. Dann ist axs keine liste sondern nur ein ding
            axs=[axs]
        #fig.subplots_adjust(hspace=0.5)
        for index,row in chunk.iterrows():
            #print(index)
            index=index*3
            reverse=False
            title="bla"
            strand=row[col_strand]
            start=int(row[col_start])-1
            stop=int(row[col_stop])-1
            total=len(row[col_ref_sequence])
            # variant position is made to be 0 based
            variantPos=int(row[col_variant_pos])-1
            if variant_pos_zero_based:
                variantPos=int(row[col_variant_pos])
                start=int(row[col_start])
                stop=int(row[col_stop])
            variantEffect = col_mpralm_logFC
            if variantEffect != "":
                variantEffect_string = f"logFC: {round(row[col_mpralm_logFC], 3)}"
                variantEffect_number = round(row[col_mpralm_logFC], 3)

            for p in range (0, len(data_base_entries), 1):
                if data_base_entries[p].split("\n")[0].split()[0] == str(row[col_motif_id]):
                    motif_name = data_base_entries[p].split("\n")[0].split()[0]
                    motif_title_string = f'Motif {motif_name}'
                    if hocomoco:
                        string_with_hocomoco_link = f'Motif <a href="https://hocomoco12.autosome.org/motif/{motif_name}#mainInfo" > <text y="{_text_size}"> Motif: {motif_name}</text> </a>'
                        motif_title_string = string_with_hocomoco_link
                    plot_logo(data_base_entries[p], axs, index+0, motif_title_string, start, stop, total, variantPos, strand, variantEffect_number, positive_variant_color, negative_variant_color, plotting_background)
                    break

            # seq 1 to pwm and plot
            ref_title = prepare_plotable_sequence_id(f"ID ({row[col_strand]}) {variantEffect_string}:\n{str(row[col_ref_id])}", row_length=30)
            plot_seq_logo(row[col_ref_sequence], axs, index+1,  ref_title, reverse, False, start, stop, variantPos, variantEffect_number, positive_variant_color, negative_variant_color, plotting_background)

            # seq 2 to pwm and plot
            alt_title = prepare_plotable_sequence_id(f"ID ({row[col_strand]}) {variantEffect_string}:\n{str(row[col_alt_id])}", row_length=30)
            plot_seq_logo(row[col_alt_sequence], axs, index+2,  alt_title, reverse, True, start, stop, variantPos, variantEffect_number, positive_variant_color, negative_variant_color, plotting_background)


        #fig.tight_layout()
        chunk_num_str = f"_{chunkNumber}" if chunkNumber > 2 else ""
        output_path = "/".join(out.split("/")[:-1])
        output_name = f'{out.split("/")[-1]}{chunk_num_str}.svg'
        plt.savefig(os.path.join(output_path, output_name))
        plt.close()

def prepare_plotable_sequence_id(string, row_length=30):
    """
    Strings reach into the sequences: only row_length characters are displayable: add "\n" after row_length characters
    """
    return "\n".join([string[i:i+row_length] for i in range(0, len(string), row_length)])


def plot_logo(data_base_entry, axs, axsrow, title, start, stop, total, variantPos, strand, variantEffect, positive_variant_color, negative_variant_color, plotting_background):
    #prepare fwd seq
    temp_PWM=[[],[],[],[]]
    for i in range (0, start, 1):
        temp_PWM[0].append(0.25)
        temp_PWM[1].append(0.25)
        temp_PWM[2].append(0.25)
        temp_PWM[3].append(0.25)

    for k in range (2, len(data_base_entry.split("\n")), 1): #2 bis -2 damit header und fuss weg fallen
        if len(data_base_entry.split("\n")[k].split()) == 4:
            temp_PWM[0].append(float(data_base_entry.split("\n")[k].split()[0]))
            temp_PWM[1].append(float(data_base_entry.split("\n")[k].split()[1]))
            temp_PWM[2].append(float(data_base_entry.split("\n")[k].split()[2]))
            temp_PWM[3].append(float(data_base_entry.split("\n")[k].split()[3]))

    for i in range (stop, total, 1):
        temp_PWM[0].append(0.25)
        temp_PWM[1].append(0.25)
        temp_PWM[2].append(0.25)
        temp_PWM[3].append(0.25)

    cwm_fwd = np.array(temp_PWM) - 0.25

    #prepare rev seq
    temp_PWM=[[],[],[],[]]
    for i in range (0, start, 1):
        temp_PWM[0].append(0.25)
        temp_PWM[1].append(0.25)
        temp_PWM[2].append(0.25)
        temp_PWM[3].append(0.25)

    for k in range (len(data_base_entry.split("\n"))-3, 1, -1): #2 bis -2 damit header und fuss weg fallen
        if len(data_base_entry.split("\n")[k].split()) == 4:
            temp_PWM[0].append(float(data_base_entry.split("\n")[k].split()[3]))
            temp_PWM[1].append(float(data_base_entry.split("\n")[k].split()[2]))
            temp_PWM[2].append(float(data_base_entry.split("\n")[k].split()[1]))
            temp_PWM[3].append(float(data_base_entry.split("\n")[k].split()[0]))

    for i in range (stop, total, 1):
        temp_PWM[0].append(0.25)
        temp_PWM[1].append(0.25)
        temp_PWM[2].append(0.25)
        temp_PWM[3].append(0.25)

    cwm_rev = np.array(temp_PWM) - 0.25

    cwm_fwd=cwm_fwd.tolist()
    cwm_rev=cwm_rev.tolist()


    df_PWM_f=pd.DataFrame()
    df_PWM_f["A"]=cwm_fwd[0]
    df_PWM_f["C"]=cwm_fwd[1]
    df_PWM_f["G"]=cwm_fwd[2]
    df_PWM_f["T"]=cwm_fwd[3]

    df_PWM_r=pd.DataFrame()
    df_PWM_r["A"]=cwm_rev[0]
    df_PWM_r["C"]=cwm_rev[1]
    df_PWM_r["G"]=cwm_rev[2]
    df_PWM_r["T"]=cwm_rev[3]

                                    #print (df_PWM)

    if strand=="+":
                                    # create forward Logo object
        PWM_logo = logomaker.Logo(df_PWM_f,
                                shade_below=.5,
                                fade_below=.5,
                                ax=axs[axsrow]
                                )

        # style using Logo methods
        PWM_logo.style_spines(visible=False)
        PWM_logo.style_spines(spines=['left'], visible=True)
        PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

        # style using Axes methods
        plot_title = title + " fwd"
        PWM_logo.ax.set_ylabel(plot_title , fontsize=3, rotation="horizontal", labelpad=10)
        PWM_logo.ax.yaxis.set_tick_params(labelsize=3)
        PWM_logo.ax.xaxis.set_ticks_position('none')
        PWM_logo.highlight_position_range(pmin=start, pmax=stop, color='silver')
        PWM_logo.highlight_position_range(pmin=variantPos, pmax=variantPos, color=positive_variant_color)
        if col_mpralm_logFC != "":
            # if variant effect positive: highlight positive_variant_color (default: cyan) else negative_variant_color (default magenta)
            if variantEffect < 0:
                PWM_logo.highlight_position_range(pmin=variantPos, pmax=variantPos, color=negative_variant_color)

        PWM_logo.ax.set_title(plot_title, fontsize=3)
        #PWM_logo.ax.xaxis.set_tick_params(pad=-1)
    elif strand=="-":
        # create rev Logo object
        PWM_logo = logomaker.Logo(df_PWM_r,
                                shade_below=.5,
                                fade_below=.5,
                                ax=axs[axsrow]
                                )

        # style using Logo methods
        PWM_logo.style_spines(visible=False)
        PWM_logo.style_spines(spines=['left'], visible=True)
        PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

        # style using Axes methods
        PWM_logo.ax.set_ylabel(title + " rev" , fontsize=3, rotation="horizontal", labelpad=10)
        PWM_logo.ax.yaxis.set_tick_params(labelsize=3)
        PWM_logo.ax.xaxis.set_ticks_position('none')
        PWM_logo.highlight_position_range(pmin=start, pmax=stop, color='silver')
        PWM_logo.highlight_position_range(pmin=variantPos, pmax=variantPos, color=positive_variant_color)
        if col_mpralm_logFC != "":
            # if variant effect positive: highlight positive_variant_color (default: cyan) else negative_variant_color (default magenta)
            if variantEffect < 0:
                PWM_logo.highlight_position_range(pmin=variantPos, pmax=variantPos, color=negative_variant_color)
        #PWM_logo.ax.xaxis.set_tick_params(pad=-1, labelsize=3)
        #PWM_logo.ax.set_title(title + " rev", fontsize=3)

def plot_seq_logo(seq, axs, axsrow, title, reverse, last, start, stop, variantPos, variantEffect, positive_variant_color, negative_variant_color, plotting_background):
    seq1PWM=[[],[],[],[]]
    letterSize=0.5
    for base in seq:
        if base == "A":
            seq1PWM[0].append(letterSize)
            seq1PWM[1].append(0)
            seq1PWM[2].append(0)
            seq1PWM[3].append(0)
        if base == "C":
            seq1PWM[0].append(0)
            seq1PWM[1].append(letterSize)
            seq1PWM[2].append(0)
            seq1PWM[3].append(0)
        if base == "G":
            seq1PWM[0].append(0)
            seq1PWM[1].append(0)
            seq1PWM[2].append(letterSize)
            seq1PWM[3].append(0)
        if base == "T":
            seq1PWM[0].append(0)
            seq1PWM[1].append(0)
            seq1PWM[2].append(0)
            seq1PWM[3].append(letterSize)


    # background = 0 #0.25
    cwm_fwd = np.array(seq1PWM)
    cwm_rev = cwm_fwd[::-1, ::-1]

    cwm_fwd=cwm_fwd.tolist()
    cwm_rev=cwm_rev.tolist()


    df_PWM_f=pd.DataFrame()
    df_PWM_f["A"]=cwm_fwd[0]
    df_PWM_f["C"]=cwm_fwd[1]
    df_PWM_f["G"]=cwm_fwd[2]
    df_PWM_f["T"]=cwm_fwd[3]

    df_PWM_r=pd.DataFrame()
    df_PWM_r["A"]=cwm_rev[0]
    df_PWM_r["C"]=cwm_rev[1]
    df_PWM_r["G"]=cwm_rev[2]
    df_PWM_r["T"]=cwm_rev[3]


    PWM_logo = logomaker.Logo(df_PWM_f,
                            shade_below=.5,
                            fade_below=.5,
                            ax=axs[axsrow],
                            #color_scheme='dimgray'
                            )

    # style using Logo methods
    PWM_logo.style_spines(visible=False)
    if last == True:
        PWM_logo.style_spines(spines=['bottom'], visible=True)
        #PWM_logo.ax.set_ylabel(title, labelpad=-1, fontsize=3)
    #else:
        #PWM_logo.ax.set_ylabel(title, labelpad=-1, fontsize=3)
        #PWM_logo.style_spines(spines=['left'], visible=True)
    PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    PWM_logo.ax.set_ylabel(title,  fontsize=3, rotation="horizontal", labelpad=20)
    # style using Axes methods

    PWM_logo.ax.yaxis.set_tick_params(left=False, labelsize=3)
    PWM_logo.ax.set(yticklabels=[])
    PWM_logo.ax.xaxis.set_ticks_position('none')
    PWM_logo.ax.xaxis.set_tick_params(pad=-1, labelsize=2)
    #PWM_logo.ax.set_title(title + " fwd", fontsize=3)
    PWM_logo.highlight_position_range(pmin=start, pmax=stop, color='silver')
    PWM_logo.highlight_position_range(pmin=variantPos, pmax=variantPos, color=positive_variant_color)
    if col_mpralm_logFC != "":
        # if variant effect positive: highlight positive_variant_color (default: cyan) else negative_variant_color (default magenta)
        if variantEffect < 0:
            PWM_logo.highlight_position_range(pmin=variantPos, pmax=variantPos, color=negative_variant_color)

    # create rev Logo object
    if reverse==True:
        PWM_logo = logomaker.Logo(df_PWM_r,
                                shade_below=.5,
                                fade_below=.5,
                                ax=axs[axsrow+1]
                                )

        # style using Logo methods
        PWM_logo.style_spines(visible=False)
        PWM_logo.style_spines(spines=['left', 'bottom'], visible=True)
        PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

        # style using Axes methods
        #PWM_logo.ax.set_ylabel("position frequency", labelpad=-1, fontsize=3)
        #PWM_logo.ax.yaxis.set_tick_params(labelsize=3)
        PWM_logo.ax.xaxis.set_ticks_position('none')
        PWM_logo.ax.xaxis.set_tick_params(pad=-1)
        PWM_logo.ax.set_title(title + " rev", fontsize=3)


if __name__ == "__main__":
    #aligned_motif_file, PWM_file, out = "MPRAlm_Alec_diffTeloHEAC_CTRL_vs_6h_withMotifs.csv", "TeloHAEC_CTRL_sigPWMs.txt", "bla"
    cli()