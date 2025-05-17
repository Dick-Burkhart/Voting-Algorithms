  
 Module Constant
   Use Precisn
   Real(Dblp), Parameter :: PI= 3.14159265358979,  Two_PI= 2.0*PI, PI_dv2= PI / 2.0
   Real(Dblp), Parameter :: deg_to_rad= PI / 180.0,  rad_to_deg= 180.0 / PI
   Real(Dblp), Parameter :: nmi_to_m  = 1852.0,    m_to_nmi= 1.0 / 1852.0 
   Real(Dblp), Parameter :: nmi_to_ft = 6076.115, ft_to_nmi= 1.0 / 6076.115
   Real(Dblp), Parameter :: Gs_to_nmisec2= 0.0053

   Real(Dblp), Parameter :: earth_radius= 6378137.0 * m_to_nmi     ! equatorial radius - nmi (WGS-84 ) 
   Real(Dblp), Parameter :: earth_rad43 = (4.0/3.0) * earth_radius ! nmi
   Real(Dblp), Parameter :: speed_of_light= 2.997925E8 * m_to_nmi       ! nmi / sec
   Real(Dblp), Parameter :: e_sq= 0.00669437999014131699613723354  ! Earth ellipticity squared (WGS-84 )
 End Module Constant
