import os
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.special import betaln
from scipy.stats import beta
from scipy.stats import binom_test
import seaborn as sns
import yaml
import family_qc
import marker_qc
import trait_predict
import generate_report


def get_args():
    from argparse import RawTextHelpFormatter
    global args
    #### Checking requirememt ##############

    examples = "python mas_pipeline.py -c [yaml_configure_file] -o [outdir, default is working dir] -s [012] \n"
    
    parse = argparse.ArgumentParser(description=f"Step0 QC of samples and markers:  -s 0  \n\
    automatically identify the individual with high missing rate and marker not in ld  \n\
    Sample needs to be removed can be manually adjusted by setting RM_IDV in the configure yaml file\n\
    Markers used for prediction can be manually adjusted under 03.marker_QC\n\
Step1 resume run after manually adjusted configures :  -s 12 \n\
    After qc, the pipeline can be resumed by assigning -o outdir, the same dir from -s 0   \n\
    After manually adjusting sample or markers run -s 12,  \n\
    do not do step 0, which is the qc again, it will DESTROY your modification \n\
Step2 is generating report only -s 2\n ", epilog=examples,formatter_class=RawTextHelpFormatter)
     
 
    parse.add_argument('-c', dest='CONFIG', type=str, required=True, help='YAML format configure file for each family')
    parse.add_argument('-o', dest='OUTDIR', type=str, required=False, help='set the PATH of the output, default is the working dir')
    parse.add_argument('-s', dest='STEPS', type=str, required=False, help='0: family and marker QC only, 1: prediction with key file, 2: generate report')


    args = parse.parse_args()
    
    if args.STEPS is None:
        args.STEPS= list("012")
    else:
        args.STEPS= list(args.STEPS)

    return args

def get_config():
    global cfg
    with open(args.CONFIG) as file:
        para_list = yaml.load(file, Loader=yaml.FullLoader)
    print(para_list)
    cfg=argparse.Namespace(**para_list)

    print(cfg.FAMILY_PATH)


def main():
    print(cfg.FAMILY_PATH)
    print(args.OUTDIR)
    if args.OUTDIR is None:
        start_time=datetime.now()
        date_time = start_time.strftime("%m%d%Y%H%M%S")
        args.OUTDIR = os.path.join(os.getcwd(), "test"+date_time)
        if "0" in args.STEPS:
            try:
                os.makedirs(args.OUTDIR)
            except FileExistsError:
                 #directory already exists
        
                print(f"{args.OUTDIR} is exist, result will be rewrite, please use a new " )
                exit(0)

    print(args)

    for dir in ["01.GT_files","02.family_QC","03.marker_QC","04.final_out","05.report"]:
        print(os.path.join(args.OUTDIR,dir))
        try:
            os.makedirs(os.path.join(args.OUTDIR,dir))
        except FileExistsError:
            pass
    print("0" in args.STEPS)
    if "0" in args.STEPS:
        os.system(f"cp {cfg.FAMILY_PATH}/hap_genotype {args.OUTDIR}/01.GT_files")
        family_qc.convert2vcf( os.path.join(args.OUTDIR,"01.GT_files","hap_genotype"),os.path.join(args.OUTDIR,"01.GT_files","hap_genotype") )
        family_qc.vcf_missindv(os.path.join(args.OUTDIR,"01.GT_files","hap_genotype.vcf"),os.path.join(args.OUTDIR,"02.family_QC","hap_genotype.vcf"))
        missing_dict=family_qc.plot_miss_hist(os.path.join(args.OUTDIR,"02.family_QC","hap_genotype.vcf.imiss"),os.path.join(args.OUTDIR,"02.family_QC","hap_genotype.vcf"),float(cfg.IDV_MISS))
        mds_dict=family_qc.vcf_mdsindv(os.path.join(args.OUTDIR,"01.GT_files","hap_genotype.vcf") ,os.path.join(args.OUTDIR,"02.family_QC","hap_genotype.vcf"))

        ### remove individual considering Missing and MDS
        oup=open(os.path.join(args.OUTDIR,"01.GT_files", "rm_missing_IDV"),"w")
        for key in missing_dict:
            if missing_dict[key] > float(cfg.IDV_MISS):
                oup.write(f"{key}\n")
                #print(f"rmmmmmmmmmmmmmmmmmmmmmmm    {key}, { missing_dict[key]}")
        oup.close()

        print("cfg.RM_IDV",cfg.RM_IDV)
        tf=os.path.join(args.OUTDIR,"01.GT_files","RM_IDV_total")
        open(tf, 'w').close()
        if cfg.RM_IDV:
            os.system("cp %s %s" % ( cfg.RM_IDV ,  os.path.join(args.OUTDIR,"01.GT_files","RM_IDV_total" )))
            os.system("cat  %s  >> %s"  %  (  os.path.join(args.OUTDIR,"01.GT_files", "rm_missing_IDV")  , os.path.join(args.OUTDIR,"01.GT_files","RM_IDV_total" ) ))
        else:
            os.system("cat  %s  >> %s"  %  (  os.path.join(args.OUTDIR,"01.GT_files", "rm_missing_IDV")  , os.path.join(args.OUTDIR,"01.GT_files","RM_IDV_total" ) ))
        cmd = "vcftools --vcf   %s   --remove %s  --recode --out %s" %  ( os.path.join(args.OUTDIR,"01.GT_files","hap_genotype.vcf")
    , os.path.join(args.OUTDIR,"01.GT_files","RM_IDV_total" ),  os.path.join(args.OUTDIR,"01.GT_files","IDV_fil") )
        print(cmd)

        ps = subprocess.Popen(cmd, shell=True)
        output = ps.communicate()[0]
        total_remove=file_len(os.path.join(args.OUTDIR,"01.GT_files","RM_IDV_total" ))
        print("total_remove",total_remove)

        total_indv= file_len(os.path.join(args.OUTDIR,"02.family_QC", "hap_genotype.vcf.imiss")) -1
        print("total_indv",total_indv)

        if (total_indv - total_remove) <2:
            print("less than two individual left")
            exit()

        family_qc.vcf_missindv(os.path.join(args.OUTDIR,"01.GT_files","IDV_fil.recode.vcf"),os.path.join(args.OUTDIR,"02.family_QC","IDV_fil.recode.vcf"))
        missing_dict2=family_qc.plot_miss_hist(os.path.join(args.OUTDIR,"02.family_QC","IDV_fil.recode.vcf.imiss"),os.path.join(args.OUTDIR,"02.family_QC","IDV_fil.recode.vcf"), float(cfg.IDV_MISS) )
        mds_dict2=family_qc.vcf_mdsindv(os.path.join(args.OUTDIR,"01.GT_files","IDV_fil.recode.vcf") ,os.path.join(args.OUTDIR,"02.family_QC","IDV_fil.recode.vcf"))

        if os.path.exists(os.path.join(args.OUTDIR, "01.GT_files", "IDV_fil.recode.vcf")):
            gt_inp = os.path.join(args.OUTDIR, "01.GT_files", "IDV_fil.recode.vcf")
        else:
            gt_inp = os.path.join(args.OUTDIR, "01.GT_files", "hap_genotype.vcf")

        for marker_file in cfg.MARKER_TRAIT_FILE:
            key_file=os.path.join(cfg.MARKER_TRAIT_PATH, marker_file)

            oup=os.path.join(args.OUTDIR,"03.marker_QC",marker_file)
            print(gt_inp,key_file,oup)
            r2_data_frame=marker_qc.marker_ld(gt_inp, key_file, oup )
            if r2_data_frame.shape[0]>1:
                marker_qc.ld_plot(r2_data_frame,oup)

            ### marker filterint  r2>0.5
            ### get rid of diagnal , if the row max r2 >0.5, keep in the new df, if no marker left, random pick one marker
            xf= r2_data_frame
            xf_masked=xf.mask(np.equal(*np.indices(xf.shape)))
            maxValuesObj = xf_masked.max(axis=1)
            print(xf_masked.index.values[maxValuesObj>0.5])
            picked_marker= xf_masked.index.values[maxValuesObj>0.5]
            if len(picked_marker) ==0:
                picked_marker=  xf_masked.index.values[list(np.random.choice(len(maxValuesObj), 1))]

            ### save the update version of key files and rerun the ploting
            fname=os.path.join(args.OUTDIR,"03.marker_QC",marker_file)
            outf=open(fname,"w")
            with open(key_file,"r") as f:
                for line in f:
                    L=line.strip().split()
                    if L[0] in picked_marker:
                        outf.write(line)
                    else:
                        outf.write("#"+line)
            outf.close()

            ### replot corrlation
            oup=os.path.join(args.OUTDIR,"03.marker_QC","filt."+marker_file)
            key_file=os.path.join(args.OUTDIR,"03.marker_QC",marker_file)
            r2_data_frame=marker_qc.marker_ld(gt_inp, key_file, oup )
            if r2_data_frame.shape[0]>1:
                marker_qc.ld_plot(r2_data_frame,oup)

            ### gt in excel  with all markers for, this step is for test only.
            ###  remove the undesired sample
            rm_ind_dic={}
            with open(os.path.join(args.OUTDIR,"01.GT_files","RM_IDV_total"),"r") as f:
                for line in f:
                    L=line.strip().split()
                    rm_ind_dic[L[0]]=1


            gt_hap=os.path.join(args.OUTDIR,"01.GT_files","hap_genotype")
            gt_hap_fil=open(os.path.join(args.OUTDIR,"01.GT_files","hap_genotype.fil"),"w")
            with open(gt_hap,"r") as f:
                for line in f:
                    L=line.strip().split()
                    #print(L)
                    if line.startswith("Locus") or line.startswith("\t"):
                        #print(L)
                        total_col=len(L)
                        pick_col=[counter for counter, value in enumerate(L) if value not in rm_ind_dic ]
                        print(pick_col)
                        gt_hap_fil.write("%s\n" % ("\t".join([ L[i]  for i in pick_col])))    
                        
                    elif len(L)== total_col:### in some row, due to too many missing data, col 2 sometimes is missing
                        gt_hap_fil.write("%s\n" % ("\t".join([ L[i]  for i in pick_col])))
            gt_hap_fil.close()


            
        for marker_file in cfg.MARKER_TRAIT_FILE:
            
            hap_inp=os.path.join(args.OUTDIR,"01.GT_files","hap_genotype.fil")
            key_file=os.path.join(cfg.MARKER_TRAIT_PATH, marker_file)
            oup=os.path.join(args.OUTDIR,"04.final_out", "gt."+marker_file+".xlsx")
            trait_predict.hap_heatmap( hap_inp, key_file , oup )


    ### plot allele frequency for all the markers( desirable, undesirable, missing)
        for marker_file in cfg.MARKER_TRAIT_FILE:
            hap_inp=os.path.join(args.OUTDIR,"01.GT_files","hap_genotype.fil")
            key_file=os.path.join(cfg.MARKER_TRAIT_PATH, marker_file)
            oup=os.path.join(args.OUTDIR,"03.marker_QC", marker_file+".freq")
            #marker_qc.allele_freq(hap_inp, key_file , oup)
            marker_qc.gt_freq(hap_inp, key_file , oup)
          


    print("1" in args.STEPS)
    if "1" in args.STEPS:
        ## define gt file
        
        ## replot ld plot according to new keys have been chosen, using vcf with undesired individual removed
        if os.path.exists(os.path.join(args.OUTDIR, "01.GT_files", "IDV_fil.recode.vcf")):
            gt_inp = os.path.join(args.OUTDIR, "01.GT_files", "IDV_fil.recode.vcf")
        else:
            gt_inp = os.path.join(args.OUTDIR, "01.GT_files", "hap_genotype.vcf")
       
 
        for marker_file in cfg.MARKER_TRAIT_FILE:
            oup=os.path.join(args.OUTDIR,"03.marker_QC","filt."+marker_file)
            key_file=os.path.join(args.OUTDIR,"03.marker_QC",marker_file)
            r2_data_frame=marker_qc.marker_ld(gt_inp, key_file, oup )
            if r2_data_frame.shape[0]>1:
                marker_qc.ld_plot(r2_data_frame,oup)

        """ 
        ## using all individual
        if os.path.exists(os.path.join(args.OUTDIR, "01.GT_files", "hap_genotype.fil")):
            hap_inp = os.path.join(args.OUTDIR, "01.GT_files", "hap_genotype.fil")
        else:
            hap_inp = os.path.join(args.OUTDIR, "01.GT_files", "hap_genotype")
        """ 

        hap_inp = os.path.join(args.OUTDIR, "01.GT_files", "hap_genotype")
        for marker_file in cfg.MARKER_TRAIT_FILE:
            print("start running prediction for", marker_file)
            key_file=os.path.join(args.OUTDIR,"03.marker_QC", marker_file)
            oup=os.path.join(args.OUTDIR,"04.final_out", "pred."+marker_file)
            bf_list=trait_predict.bf_hap_predict(hap_inp, key_file , oup )
            trait_predict.plot_bf_hist( bf_list , oup+".bf_hist.png" )


    if "2" in args.STEPS:
        print("-----------------start generating report--------------\n")
        ouf=open(os.path.join(args.OUTDIR, "05.report", "multiqc_config.yaml"),"w")
        ouf.write("for data_format: %r\n\n" % 'yaml')
    

        ouf.write("title: %r\n" % "Marker assisted selection generated by Vitisgen2")
        ouf.write("subtitle: %r\n" % "a family  with %s individual were included and %s markers-trait were predicted" % ( 
           file_len(os.path.join(args.OUTDIR,"02.family_QC", "hap_genotype.vcf.imiss")) -1 ,
          len(cfg.MARKER_TRAIT_FILE)   ))
        #"intro: "You can find out more about MultiQC at <a href="http://multiqc.info">multiqc.info</a>"
        #custom_message: "Remember - MultiQC is awesome!"
        ouf.write("report_header_info:\n")
        #ouf.write("    - Contact E-mail: %r\n" % 'phil.ewels@scilifelab.se')
        ouf.write("    - Application Type: %r\n" % 'rhAmpSeq')
        ouf.write("    - Sequencing Platform: %r\n" %'NextSeq500 ')
        ouf.write("    - Sequencing Setup: %r\n" % '2x150bp')
        ouf.write("custom_logo:  %r\n" % "./vitisgen.png")

        ouf.write("table_cond_formatting_rules: \n") 
        for column in [f"pred.{marker_file}.bf.txt" for marker_file in cfg.MARKER_TRAIT_FILE]:
            ouf.write("            %s:\n" % column.replace(".","_"))
            ouf.write("                pass:\n")
            ouf.write("                    - gt: 0\n")
            ouf.write("                warn:\n")
            ouf.write("                    - eq: 0\n")
            ouf.write("                fail:\n")
            ouf.write("                    - lt: 0\n")

        column="IDVrm"
        if 1:
            ouf.write("            %s:\n" % column.replace(".","_"))
            ouf.write("                pass:\n")
            ouf.write("                    - gt: 0\n")
            ouf.write("                warn:\n")
            ouf.write("                    - eq: 0\n")
            ouf.write("                fail:\n")
            ouf.write("                    - lt: 0\n")


        ouf.write("custom_content:\n")
        ouf.write("    order:\n")

      

        # 
        #Step 1 move result figures and generate yaml tables for family QC
        # 
        #1.1 get individual missing qc 
        custom_data_list=[] ## generate custom_data_lines

        sp_list=[] ## genetate custom sp lines

        missing_plot1=os.path.join(args.OUTDIR, "02.family_QC", "hap_genotype.vcfmissing_indv.png")
        if os.path.exists(missing_plot1) :
            os.system("cp %s %s" % (missing_plot1, os.path.join(args.OUTDIR, "05.report", "missing_before_filter_mqc.png")))
            ouf.write("        - missing_before_filter\n")
            custom_data_lines=[]
            custom_data_lines.append("    missing_before_filter:\n")
            custom_data_lines.append("        id: %r\n" %  'missing_before_filter')
            custom_data_lines.append("        section_anchor:  %r\n" %  'missing_before_filter')
            custom_data_lines.append("        section_name: %r\n" % 'missing before filter')
            custom_data_lines.append("        description : %r\n" % 'missing rate for each sample before filter')
            custom_data_list.append(custom_data_lines)             

            sp_line=[]
            sp_line.append("    missing_before_filter:\n")
            sp_line.append("        fn: %r\n" %  'missing_before_filter_mqc.png')
            sp_list.append(sp_line)

        missing_plot2=os.path.join(args.OUTDIR, "02.family_QC", "IDV_fil.recode.vcfmissing_indv.png")
        if os.path.exists(missing_plot2):
            os.system("cp %s %s" % (missing_plot2, os.path.join(args.OUTDIR, "05.report", "missing_after_filter_mqc.png")))
            ouf.write("        - missing_after_filter\n")
            custom_data_lines=[]
            custom_data_lines.append("    missing_after_filter:\n")
            custom_data_lines.append("        id: %r\n" %  'missing_after_filter')
            custom_data_lines.append("        section_anchor: %r\n" %  'missing_after_filter')

            custom_data_lines.append("        section_name: %r\n" % 'missing after filter')
            custom_data_lines.append("        description : %r\n" % 'missing rate for each sample after filter')
            custom_data_list.append(custom_data_lines)             

            sp_line=[]
            sp_line.append("    missing_after_filter:\n")
            sp_line.append("        fn: %r\n" %  'missing_after_filter_mqc.png')
            sp_list.append(sp_line)

  
        #1.2 get individual mds 
        # from 01.GT_files get error of sample names
        mds_out1=os.path.join(args.OUTDIR, "02.family_QC", "hap_genotype.vcf.mds")
        if  os.path.exists(mds_out1):
            oup=os.path.join(args.OUTDIR, "05.report", "mds_before_filter_mqc.yaml")
            para_dic={"id":"mds_before_filter",
                    "name":"mds_before_filter",
                    "description":"MDS calculated by plink",
                    "plot_type":"scatter",
                    "pconfig_id":"mqc_mds_scatter_plot",
                    "pconfig_title":"MDS plot",
                    "pconfig_x":"Component 1",
                    "pconfig_y":"Component 2"
                     }
           
            generate_report.pca_yaml(mds_out1,para_dic,cfg.MOTHER_NAMES ,cfg.FATHER_NAMES,oup)
            ouf.write("        - mds_before_filter\n")
            custom_data_lines=[]
            custom_data_lines.append("    mds_before_filter:\n")
            custom_data_lines.append("        id: %r\n" %  'mds_before_filter')
            custom_data_lines.append("        section_anchor: %r\n" %  'mds_before_filter')
            custom_data_lines.append("        section_name: %r\n" % 'mds_before_filter')
            custom_data_lines.append("        description : %r\n" % 'MDS plot for each sample before filter,red for mother, blue for father and purple for offsprings')
            custom_data_list.append(custom_data_lines)             

            sp_line=[]
            sp_line.append("    mds_before_filter:\n")
            sp_line.append("        fn: %r\n" %  'mds_before_filter_mqc.yaml')
            sp_list.append(sp_line)

       
        mds_out2=os.path.join(args.OUTDIR, "02.family_QC", "IDV_fil.recode.vcf.mds")
        if  os.path.exists(mds_out1):
            oup=os.path.join(args.OUTDIR, "05.report", "mds_after_filter_mqc.yaml")
            para_dic={"id":"mds_after_filter",
                    "name":"mds_after_filter",
                    "description":"MDS calculated by plink",
                    "plot_type":"scatter",
                    "pconfig_id":"mqc_mds_scatter_plot",
                    "pconfig_title":"MDS plot",
                    "pconfig_x":"Component 1",
                    "pconfig_y":"Component 2"
                     }
           
            generate_report.pca_yaml(mds_out2,para_dic,cfg.MOTHER_NAMES ,cfg.FATHER_NAMES,oup)
            ouf.write("        - mds_after_filter\n")
            custom_data_lines=[]
            custom_data_lines.append("    mds_after_filter:\n")
            custom_data_lines.append("        id: %r\n" %  'mds_after_filter')
            custom_data_lines.append("        section_anchor: %r\n" %  'mds_after_filter')
            custom_data_lines.append("        section_name: %r\n" % 'mds_after_filter')
            custom_data_lines.append("        description : %r\n" % 'MDS plot for each sample after filter, red for mother, blue for father and purple for offsprings')
            custom_data_list.append(custom_data_lines)             

            sp_line=[]
            sp_line.append("    mds_before_filter:\n")
            sp_line.append("        fn: %r\n" %  'mds_after_filter_mqc.yaml')
            sp_list.append(sp_line)
            
        
        
        for marker_file in cfg.MARKER_TRAIT_FILE:
            
            corr_plot1=os.path.join(args.OUTDIR, "03.marker_QC", f"{marker_file}.r2.png")
            if os.path.exists(corr_plot1) :
                os.system("cp %s %s" % (corr_plot1, os.path.join(args.OUTDIR, "05.report", f"{marker_file}_corr_before_mqc.png")))
               
                replace_id=marker_file.replace(".","_")
                ouf.write("        - %s\n" % (f"{replace_id}_corr_before"  ))
                custom_data_lines=[]
                custom_data_lines.append("    %s:\n" % ( f"{replace_id}_corr_before" ))
                custom_data_lines.append("        id: %r\n" % ( '%s' %  f"{replace_id}_corr_before"))
                custom_data_lines.append("        section_anchor: %r\n" % ( '%s' %   f"{replace_id}_corr_before"))
                custom_data_lines.append("        section_name: %r\n" %  f"{replace_id}_corr_before")
                custom_data_lines.append("        description : %r\n" % 'correlation plot for each marker before filtering')
                custom_data_list.append(custom_data_lines) 
                sp_line=[]
                sp_line.append("    %s:\n" % ( f"{replace_id}_corr_before" ))
                sp_line.append("        fn: %r\n" % ('%s' % f"{marker_file}_corr_before_mqc.png"))
                sp_list.append(sp_line)

            frq_plot=os.path.join(args.OUTDIR, "03.marker_QC", f"{marker_file}.freq")
            if os.path.exists(frq_plot) :
                os.system("cp %s %s" % (frq_plot, os.path.join(args.OUTDIR, "05.report", f"{marker_file}_freq_mqc.txt")))
                replace_id=marker_file.replace(".","_")
                ouf.write("        - %s\n" % (f"{replace_id}_freq"  ))
                custom_data_lines=[]
                custom_data_lines.append("    %s:\n" % ( f"{replace_id}_freq" ))
                custom_data_lines.append("        id: %r\n" % ( '%s' %  f"{replace_id}_freq"))
                custom_data_lines.append("        section_anchor: %r\n" % ( '%s' %   f"{replace_id}_freq"))
                custom_data_lines.append("        section_name: %r\n" %  f"{replace_id}_freq")
                custom_data_lines.append("        description : %r\n" % 'frequency plot for each marker before filtering')
                custom_data_list.append(custom_data_lines) 
                sp_line=[]
                sp_line.append("    %s:\n" % ( f"{replace_id}_freq" ))
                sp_line.append("        fn: %r\n" % ('%s' % f"{marker_file}_freq_mqc.txt"))
                sp_list.append(sp_line)
     


 
            corr_plot2=os.path.join(args.OUTDIR, "03.marker_QC", f"filt.{marker_file}.r2.png")
            if os.path.exists(corr_plot2):
                os.system("cp %s %s" % (corr_plot2, os.path.join(args.OUTDIR, "05.report", f"{marker_file}_corr_after_mqc.png")))
                replace_id=marker_file.replace(".","_")
                ouf.write("        - %s\n" % ( f"{replace_id}_corr_after" ))
                custom_data_lines=[]
                custom_data_lines.append("    %s:\n" % ( f"{replace_id}_corr_after" ))
                custom_data_lines.append("        id: %r\n" % ( '%s' %  f"{replace_id}_corr_after"))
                custom_data_lines.append("        section_anchor: %r\n" % ( '%s'%  f"{replace_id}_corr_after"))
                custom_data_lines.append("        section_name: %r\n" % f"{replace_id}_corr_after")
                custom_data_lines.append("        description : %r\n" % 'correlation plot for each marker after filtering')
                custom_data_list.append(custom_data_lines) 
                sp_line=[]
                sp_line.append("    %s:\n" % ( f"{replace_id}_corr_after" ))
                sp_line.append("        fn: %r\n" % ('%s' % f"{marker_file}_corr_after_mqc.png"))
                sp_list.append(sp_line)

        
            bf_hist=os.path.join(args.OUTDIR, "04.final_out", f"pred.{marker_file}.bf_hist.png")
            if os.path.exists(bf_hist) :
                os.system("cp %s %s" % (bf_hist, os.path.join(args.OUTDIR, "05.report", f"{marker_file}_bf_hist_mqc.png")))
                replace_id=marker_file.replace(".","_")
                ouf.write("        - %s\n" % ( f"{replace_id}_bf_hist" ))
                custom_data_lines=[]
                custom_data_lines.append("    %s:\n" % ( f"{replace_id}_bf_hist" ))
                custom_data_lines.append("        id: %r\n" % ( '%s' %  f"{replace_id}_bf_hist"))
                custom_data_lines.append("        section_anchor: %r\n" % ( '%s' %  f"{replace_id}_bf_hist"))
                custom_data_lines.append("        section_name: %r\n" % f"{replace_id}_bf_hist")
                custom_data_lines.append("        description : %r\n" % 'histgram of log10 BF' )
                custom_data_list.append(custom_data_lines) 
                sp_line=[]
                sp_line.append("    %s:\n" % ( f"{replace_id}_corr_before" ))
                sp_line.append("        fn: %r\n" % ('%s' % f"{marker_file}_bf_hist_mqc.png"))
                sp_list.append(sp_line)


        para_dic={"id": "Prediction_Summary",
                    "name":"Prediction_Summary",
                    "description":"summary of missing and BF of prediction",
                    "plot_type":"table"
            }

        bf_file_list=[ os.path.join(args.OUTDIR, "04.final_out", f"pred.{marker_file}.bf.txt") for marker_file in cfg.MARKER_TRAIT_FILE ]
        missing_file=os.path.join(args.OUTDIR,"02.family_QC", "IDV_fil.recode.vcf.imiss")
        oup=os.path.join(args.OUTDIR, "05.report", "all_sum_table_mqc.txt")
        generate_report.sum_table_txt(para_dic,oup,missing_file,bf_file_list,cfg.RM_IDV)
        ouf.write("        - Prediction_Summary\n")
        custom_data_lines=[]
        custom_data_lines.append("    Prediction_Summary:\n")
        custom_data_lines.append("        id: %r\n" %  'Prediction_Summary')
        custom_data_lines.append("        section_anchor: %r\n" %  'Prediction_Summary')
        custom_data_lines.append("        section_name: %r\n" % 'Prediction_Summary')
        custom_data_lines.append("        description : %r\n" % 'Summary table for each sample with missing rate and BF prediction for each trait.')
        custom_data_lines.append("        format: %r\n" % 'tsv')
        custom_data_lines.append("        plot_type: %r\n" % 'table')
       
        
        custom_data_list.append(custom_data_lines)             

        sp_line=[]
        sp_line.append("    Prediction_Summary:\n")
        sp_line.append("        fn: %r\n" %  'all_sum_table_mqc.txt')
        sp_list.append(sp_line)

        ouf.write("custom_data:\n")
        for i in custom_data_list:
            for j in i:
                ouf.write(j)
       
        ouf.write("sp:\n")
        for i in sp_list:
            for j in i:
                ouf.write(j)

          
        generate_report.sum_xlsx(args,cfg)


def file_len(file):
    return len(open(file).readlines())


if __name__ == "__main__":
    get_args()
    get_config()
    main()


