

import numpy




def calculate_consensus_genotypes(allele_counts_matrix,lower_threshold=0.2,upper_threshold=0.8):
    
    num_sites, num_samples, num_alleles = allele_counts_matrix.shape
    
    depths = allele_counts_matrix.sum(axis=2) # 2D
    freqs = allele_counts_matrix[:,:,0]*1.0/(depths+(depths==0)) # 2D
    passed_sites_matrix = (depths>0)*numpy.logical_or(freqs<=lower_threshold,freqs>=upper_threshold)
    # consensus approximation
    genotype_matrix = numpy.around(freqs)*passed_sites_matrix
    
    # genotype_matrix is a 2D data object
    
    return genotype_matrix, passed_sites_matrix

#####################################################################


#####################################################################
def calculate_unbiased_sigmasquared(allele_counts_1, allele_counts_2):
    
    '''
    Adapted from Garud, Good et al.
    Line 202 in: https://github.com/benjaminhgood/microbiome_evolution/blob/4015ec03e2cbdbfb7d062415fb6a8c869f37a32e/diversity_utils.py
    '''
    # An alternate version of a standard measure of linkage disequilibrium:
    #
    # sigma_squared= E[X]/E[Y], where X=(p_ab-pa*pb)^2 and Y=(pa*(1-pa)*pb*(1-pb))
    # rsquared=E[X/Y]
    # see McVean 2002 for more notes on the difference. 
    #
    # where we have corrected for finite sample effects

    #genotypes_1, passed_sites_1 = calculate_consensus_genotypes(allele_counts_1)
    #genotypes_2, passed_sites_2 = calculate_consensus_genotypes(allele_counts_2)
    
    
    # this asks which pairs of sites have depths >0 at BOTH sites
    # None here takes the product of the elements in the two vectors and returns a matrix. 
    joint_passed_sites = (passed_sites_1)[None,:,:]*(passed_sites_2)[:,None,:]
    # sites x sites x samples matrix
    
    # allele counts
    ns = joint_passed_sites.sum(axis=2)
    
    n11s = ((genotypes_1[None,:,:])*(genotypes_2[:,None,:])*joint_passed_sites).sum(axis=2)
    n10s = (genotypes_1[None,:,:]*(1-genotypes_2[:,None,:])*joint_passed_sites).sum(axis=2)
    n01s = ((1-genotypes_1[None,:,:])*(genotypes_2[:,None,:])*joint_passed_sites).sum(axis=2)
    n00s = ((1-genotypes_1[None,:,:])*(1-genotypes_2[:,None,:])*joint_passed_sites).sum(axis=2)
    

    
    # First calculate numerator
    rsquared_numerators = n11s*(n11s-1)*n00s*(n00s-1)
    rsquared_numerators -= 2*n10s*n01s*n11s*n00s
    rsquared_numerators += n10s*(n10s-1)*n01s*(n01s-1)
    
    #print "Before divide:"
    #print rsquared_numerators
    
    rsquared_numerators = rsquared_numerators*(ns>3.5)*1.0/(ns*(ns-1)*(ns-2)*(ns-3)+10*(ns<3.5))
    
    #print "After divide:"
    #print rsquared_numerators
    #print "---"
    # Now calculate denominator 
    # (more annoying... there are 16 terms rather than 4, so we will write them separately)
    
    #1
    rsquared_denominators = n10s*(n10s-1)*n01s*(n01s-1)    
    #2
    rsquared_denominators += n10s*n01s*(n01s-1)*n00s        
    #3
    rsquared_denominators += n10s*(n10s-1)*n01s*n11s         
    #4
    rsquared_denominators += n10s*n01s*n11s*n00s           
    #5
    rsquared_denominators += n10s*(n10s-1)*n01s*n00s       
    #6
    rsquared_denominators += n10s*n01s*n00s*(n00s-1)
    #7
    rsquared_denominators += n10s*(n10s-1)*n11s*n00s
    #8
    rsquared_denominators += n10s*n11s*n00s*(n00s-1)
    #9
    rsquared_denominators += n10s*n01s*(n01s-1)*n11s
    #10
    rsquared_denominators += n01s*(n01s-1)*n11s*n00s
    #11
    rsquared_denominators += n10s*n01s*n11s*(n11s-1)
    #12
    rsquared_denominators += n01s*n11s*(n11s-1)*n00s
    #13
    rsquared_denominators += n10s*n01s*n11s*n00s
    #14
    rsquared_denominators += n01s*n11s*n00s*(n00s-1)
    #15
    rsquared_denominators += n10s*n11s*(n11s-1)*n00s
    #16
    rsquared_denominators += n11s*(n11s-1)*n00s*(n00s-1)
    
    # divide by sample size
    rsquared_denominators = rsquared_denominators*(ns>3.5)*1.0/(ns*(ns-1)*(ns-2)*(ns-3)+10*(ns<3.5))
    
    return rsquared_numerators, rsquared_denominators


