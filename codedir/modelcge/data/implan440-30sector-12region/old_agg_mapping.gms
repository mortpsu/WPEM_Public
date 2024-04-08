SET
* Aggregate sector subsets
        s_a           Aggregate Sectors /
                        agrs    "Agricultural sectors"
                        crps    "Crop sectors"
                        eles    "Electricity sectors"
                        mins    "Mining sectors"
                        pens    "Primary energy sectors"
                        sens    "Secondary energy sectors"
                        cons    "Construction sectors"
                        mans    "Manufacturing sectors"
                        srvs    "Services sectors"
                        trns    "Transportation sectors" /

        s_map(s_a,s)  Aggregate Sector Map  /
                        agrs.( apa, cba, frs, grn, oca, osa, pfb, vna ),
                        crps.( cba, grn, oca, osa, pfb, vna ),
                        eles.( ele ),
                        mins.( min ),
                        pens.( coa, cru ),
                        sens.( ngd, per ),
                        cons.( con ),
                        mans.( bom, cem, cpm, fbm, pfm, tec, trm, wpm ),
                        srvs.( bos, fin, hlt, pub, rtl, tel ),
                        trns.( trn ) /
;

Set
* Region subsets

        cal(r)        California                 / cal /
        rwc(r)        Rest of Western states     / arz, col, ido, mta, nmo,
                                                   nva, org, uth, was, wyo /
        rus(r)        Rest of US states          / rus /

* Aggregate region subsets
        r_a           Aggregate region /
                         CA      "California"
                         ROWECC  "Rest of Western states"
                         ROUS    "Rest of U.S. States" /

        r_map(r_a,r)  Aggregate region map  /
                        CA.( cal ),
                        ROWECC.( arz, col, ido, mta, nmo, nva, org, uth, was, wyo ),
                        ROUS.( rus ) /

        r_psm(r,r_a)  Reverse of r_map
;

* Define reverse of region map
option  r_psm < r_map
