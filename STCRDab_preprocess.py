

def stcrdab_with_vab_ag(strdab_summary_filename):
    inf=open(wdir+'\\'+strdab_summary_filename)
    outf=open(wdir+'\\'+strdab_summary_filename.replace('.tsv','.abTCRMH1Ag.tsv'),'w')
    fields=['pdb', 'Bchain', 'Achain', 'Dchain', 'Gchain', 'TCRtype', 'model', 'antigen_chain',
            'antigen_type', 'antigen_name', 'antigen_het_name', 'mhc_type', 'mhc_chain1',
            'mhc_chain2', 'docking_angle', 'beta_subgroup', 'alpha_subgroup', 'gamma_subgroup',
            'delta_subgroup', 'short_header', 'date', 'compound', 'beta_organism', 'alpha_organism',
            'gamma_organism', 'delta_organism', 'antigen_organism', 'mhc_chain1_organism',
            'mhc_chain2_organism', 'authors', 'resolution', 'method', 'r_free', 'r_factor',
            'affinity', 'affinity_method', 'affinity_temperature', 'affinity_pmid', 'engineered']
    for i,line in enumerate(inf):
        if i==0:
            outf.write(line)
            continue
        info=line.strip().split('\t')
        pdbid=info[0].strip()
        dinfo=dict(zip(fields,info))
        if not dinfo['TCRtype']=='abTCR':
            continue
        if len(dinfo['antigen_chain'])==0:
            continue
        if not dinfo['mhc_type']=='MH1':
            continue
        outf.write(line)
    inf.close()
    outf.close()

wdir='C:\\Users\\Norman.Goodacre\\3T\\STCRDab'
strdab_summary_filename='tcrdb_summary_all.tsv'
stcrdab_with_vab_ag(strdab_summary_filename)
              
    
        
