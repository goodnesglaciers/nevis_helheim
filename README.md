# nevis_helheim
The Hewitt (2013) subglacial hydrology model adapted for Helheim Glacier, East Greenland.

Citations:

Stevens, L. A., Nettles, M., Davis, J. L., Creyts, T. C., Kingslake, J., Hewitt, I. J., and A. Stubblefield (2021). High meltwater throughput limits Helheim Glacier response to lake drainage, In review.

Stevens, L. A., Hewitt, I., Das, S. B., and Behn, M. D., (2018), Relationship between Greenland Ice Sheet surface speed and modeled effective pressure, Journal of Geophysical Research: Earth Surface, Journal of Geophysical Research: Earth Surface, 123(9), p2258–2278, doi:10.1029/2017JF004581.

Hewitt I., (2013), Seasonal changes in ice sheet motion due to melt water lubrication, Earth Planet. Sci. Lett., 371, p16–25, doi:10.1016/j.epsl.2013.04.022.

Model run files are given for:
    
    (1) A 2007/1–365 run with distributed runoff input: nevis_h22222_ubspatial_R365.m.
    
    (2) An M_229 lake-drainage event: nevis_h22222_ubspatial_R227_lakerampM_4tiles_Ks100_s1e6H_repo.m
    
    (3) An M_winter lake-drainage event: nevis_h22222_ubspatial_R67_lakerampM_4tiles_Ks100_s1e6H_repo.m 
    
All model code is located within the "nevis" directory. To run the model for other sheet permeability (K_s) and englacial void fraction (sigma) parameter values, edit the run file where indicated. Additional information on running the model is provided in nevis.pdf (written by Ian Hewitt).

Contact LAS with questions (laura.stevens@earth.ox.ac.uk).
