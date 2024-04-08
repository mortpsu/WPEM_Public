set 
		econr         PSM state regions 
                 /
                        1        California
                        2        Arizona
                        3        Colorado
                        4        Idaho
                        5        Montana
                        6        Nevada
                        7        New Mexico
                        8        Oregon
                        9        Utah
                        10       Washington
                        11       Wyoming
                 /
				 
	zones         PSM zone regions
                 /
                         1        New Mexico
                         2        Arizona
                         3        Nevada
                         4        California
                         6        Utah
                         7        Oregon\Washington
                         8        Montana
                         9        Idaho
                         12       Wyoming
                         13       Colorado
                 /		 
		
		zones_map(zones,r)  Aggregate region map
                 /
                        4  . cal
                        2  . arz
                        13 . col
                        9  . ido
                        8  . mta
                        1  . nmo
                        3  . nva
                        7  . ( org, was )
                        6  . uth
                        12 . wyo
                 /
		 
		econz_map (econr,zones)  Aggregate region map
                 /
                        1  . 4
                        2  . 2
                        3  . 13
                        4  . 9
                        5  . 8
                        6  . 3
                        7  . 1
                        8  . 7
                        9  . 6
                        10 . 7
                        11 . 12
                 /				 
				 
        econr_map(econr,r)  Aggregate region map
                 /
                        1  . cal
                        2  . arz
                        3  . col
                        4  . ido
                        5  . mta
                        6  . nva
                        7  . nmo
                        8  . org
                        9  . uth
                        10 . was
                        11 . wyo
                 /

        econzr_map (econr,zones,r)  Aggregate region map
                 /
                        1  .  4  . cal
                        2  .  2  . arz
                        3  .  13 . col
                        4  .  9  . ido
                        5  .  8  . mta
                        6  .  3  . nva
                        7  .  1  . nmo
                        8  .  7  . org
                        9  .  6  . uth
                        10 .  7  . was
                        11 .  12 . wyo
                 /

        r_psm_zones(r,zones)  Reverse of zones_map
        r_psm_econr(r,econr)  Reverse of state_map
;

*                Define reverse of region map
option  r_psm_zones < zones_map
option  r_psm_econr < econr_map				 