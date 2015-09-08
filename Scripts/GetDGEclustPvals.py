import subprocess

species_list = ['KRAUS', 'MOEL', 'UNC', 'WILD']
species_starts = { 'KRAUS' : 0, 'MOEL' : 4, 'UNC' : 7, 'WILD' : 11 }
for x in range(len(species_list)-1): 
  for y in range(x+1,len(species_list)):
    #calculate pvalues, pooled by species
    subprocess.call(["pvals -t0 1000 -tend 10000 -i By_species/_clust -o %s -g1 %s -g2 %s" % ('By_species/' + '_'.join((species_list[x], species_list[y], 'pvals.txt')), x, y)], shell=True)
    
    #calculate pvalues pooled by leaf
    subprocess.call(["pvals -t0 1000 -tend 10000 -i By_leaf/_clust -o %s -g1 %s -g2 %s" % ('By_leaf/' + ''.join(('All', str(x+1), str(y+1), '_pvals.txt')), x, y)], shell=True)
    
    #calculate pvalues for individual comparisons
    for species in species_list:  
      if species == 'MOEL' and y == 3:  
        pass #there is no MOEL4 sample
      else:
        subprocess.call(["pvals -t0 1000 -tend 10000 -i All_by_All/_clust -o %s -g1 %s -g2 %s" % ('All_by_All/' + ''.join((species, str(x+1), str(y+1), '_pvals.txt')), x + species_starts[species], y + species_starts[species])], shell=True)

