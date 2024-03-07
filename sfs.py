import allel
import argparse
import numpy as np
from pathlib import Path


def main_arguments():
    parser = argparse.ArgumentParser(description="Get folded SFS from VCF")
    parser.add_argument('-v', '--vcf', help='vcf',required=True)
    parser.add_argument('-o', '--outpath', help='outpath',required=True)
    parser.add_argument('-s', '--samples', help='samples to use. otherwise all samples included.',required=False)
    parser.add_argument('--blueprint_outpath', help='if using for stairwayplot. if used, L, popid, nseq, whether_folded, project_dir, mu, year_per_generation, and plot_title are required', required=False)
    parser.add_argument('--popid', required=False)
    parser.add_argument('--nseq', required=False)
    parser.add_argument('--whether_folded', required=False)
    parser.add_argument('--project_dir', required=False)
    parser.add_argument('--mu', required=False)
    parser.add_argument('--year_per_generation', required=False)
    parser.add_argument('--plot_title', required=False)
    parser.add_argument('--L', required=False)
    args = parser.parse_args()
    return(args)


def read_genotype_array(vcf,samples):

    '''
    Read only GT data from VCF file.
    Convert to Genotype array class & return
    '''

    #samples = kwargs.get('samples', None)
    
    if samples:
        genotypes = allel.read_vcf(vcf, samples=samples, fields=['calldata/GT'])
    else:
        genotypes = allel.read_vcf(vcf, fields=['calldata/GT'])

    genotypes_field = genotypes['calldata/GT']
    genoarr = allel.GenotypeArray(genotypes_field)

    return(genoarr)


def get_sfs(genoarr):

    '''
    Calculate SFS given a genotype array.
    Folded default, optional unfolded with argument:
    unfolded=True
    '''

    #unfolded = kwargs.get('unfolded', None)

    ac = genoarr.count_alleles()

    #if unfolded:
     #   sfs=allel.sfs(ac)
    #else:
    sfs = allel.sfs_folded(ac)

    return(sfs)


def write_sfs(sfs, outpath):

    '''
    Write sfs to file of choice.
    '''

    with open(outpath,'w') as f:
        f.write(' '.join(map(str,(list(sfs)))))


def write_blueprint(sfs, outpath, L, popid, nseq, whether_folded, project_dir, mu, year_per_generation, plot_title):

    '''
    Create blueprint file for use in stairwayplot
    '''
    first=str(round(int(nseq)/2))
    second=str(round((int(nseq)-2)/4))
    third=str(round((int(nseq)-2)/2))
    fourth=str(round((int(nseq)-2)*0.75))
    fifth=str(round(int(nseq)-2))
    with open(outpath,'w') as f:
        f.write("#input setting\npopid: {}\nnseq: {}\nL: {}\nwhether_folded: {}\nSFS: {}\n#smallest_size_of_SFS_bin_used_for_estimation: 1\n#largest_size_of_SFS_bin_used_for_estimation: {}\npct_training: 0.67\nnrand: {} {} {} {}\nproject_dir: {}\nstairway_plot_dir: stairway_plot_es\nninput: 200\nrandom_seed: 672\n#output setting\nmu: {}\nyear_per_generation: {}\n#plot setting\nplot_title: {}\nxrange: 0,0\nyrange: 0,0\nxspacing: 2\nyspacing: 2\nfontsize: 12".format(popid,nseq,L,whether_folded,sfs,first,second,third,fourth,fifth, project_dir,mu,year_per_generation,plot_title))


def main():

    args = main_arguments()
    vcf = args.vcf
    outpath = args.outpath

    #create list of samples to pull from VCF if requested
    if args.samples:
        samplesstr = args.samples
        samples = samplesstr.split(',')

        #if not all(isinstance(i, str) for i in samples):
            #for index, name in enumerate(samples):
                #list[index] = str(name)

        genoarr=read_genotype_array(vcf, samples)

    else:
        #read in GT fields of VCF and return array
        genoarr=read_genotype_array(vcf, None)

    #calculate sfs from GT array
    #if unfolded:
     #   sfs=get_sfs(genoarr, unfolded=True)
    #else:
    sfs=get_sfs(genoarr)

    #write sfs to file
    write_sfs(sfs, outpath)

    if args.blueprint_outpath:
        sfs=' '.join(map(str,(list(sfs))))
        print(sfs)
        #if unfolded:
         #   whether_folded="false"
        #else:
        whether_folded="true"

        write_blueprint(sfs, args.blueprint_outpath, args.L, args.popid, args.nseq, whether_folded, args.project_dir, args.mu, args.year_per_generation, args.plot_title)


if __name__ == "__main__":
    main()
