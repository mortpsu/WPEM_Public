$title Input file for PCHES power system model (PSM) - Retirement SHOCK TrialImportData

$ontext
    This GAMS program loads the data necessary for the DCOF and Unit Commitment PSM models. 
    This file was orginally created my Vijay Kumar and Dr. Mort Webster. 
    Additional edits and annotations have been added by Joseph Perla:
    
    Updated by JMP629: 2021-03-14
$offtext

* ---------------------------------------------------------------------------- *
*       model information:
* ---------------------------------------------------------------------------- *

$ontext

Number of generators 3569

Number of buses 312

Number of transmission lines 36

state name      zones id
New Mexico         1
Arizona            2
Nevada             3
California         4
CFE                5
Utah               6
Oregon             7
Washington         7
Montana            8
Idaho              9
British Columbia   10
Alberta            11
Wyoming            12
Colorado           13

fuel type name            fuel type id     number of generators 
Coal-PC                             1         131
Wind-Onshore                        2         249
DG-BTM                              3         0
Gas CT                              5         1547
Hydro                               4         621
Gas CCGT                            6         183
Steam Turbine                       7         111
Solar PV                            8         31
	-Distributed Utility-Fixed Tilt 
Solar PV                            9         4          
	-Distributed Utility-Tracking   
DC-Intertie                         10        8
DR                                  11        0
Nuclear                             12        6
PS-Hydro                            13        55
Biogas-Other                        14        276
SolarThermal                        15        13
Petroleum                           16        47
Geothermal                          17        159
ES-Generic                          18        0
EE                                  19        0
MotorLoad                           20        128

- SREM ---- regional economic model
state name      GenArea id
California         1
Arizona            2
Colorado           2
Idaho              2
Montana            2
Nevada             2
New Mexico         2
Oregon             2
Utah               2
Washington         2
Wyoming            2
Alberta            2
British Columbia   2
CFE                2

- DREM ---- regional economic model
state name      GenArea id 
California         1
Arizona            2
Colorado           3
Idaho              4
Montana            5
Nevada             6
New Mexico         7
Oregon             8
Utah               9
Washington         10
Wyoming            11
Alberta            12
British Columbia   13
CFE                14


$offtext

* ---------------------------------------------------------------------------- *
*       define model set for the unitcommitment and DCOFP programs:
* ---------------------------------------------------------------------------- *

sets
*   main sets
       i                     bus
       g                     generator
       d                     hours                               / d1*d168 /
       line                  transmission lines                  / l1*l654 /
       
*   secondary sets
       lines                 transmission lines between zones    / 1*36 /
       ft                    fuel type
       cl                    line charactoristics (reactance and thernal limit) 
                                                                 / r, x1, lm, can, current, cost /                      
*   zone and region sets   
	   zones                  modeling zones (similar to states)
*	   GenLocation            economic modeling regions (states)  / 1*%num_state% /
*	   zonesforEcon           economic modeling regions (states)  / 1*%num_state% /
	   econr                  economic modeling regions (states)

*   crosswalk sets        
	   ig(i,g)                generator to bus crosswalk  (g in i)
	   gOff(g)                generator affected by WBM   (g in WBM)
	   Gentozone(g,zones)     generator to zone crosswalk (g in zones)
	   
	   gzones(g,zones)        generator to zones crosswalk  (g in zones)
	   geconr(g,econr)        generator to economic region crosswalk  (g in econr)
	   gi(g,i)                generator to bus crosswalk (g in i)
	   gftype(g,ft)           generator to fuel type crosswalk  (g in i)
;

* ---------------------------------------------------------------------------- *
*       define model generator parameters the unitcommitment and DCOFP programs:
* ---------------------------------------------------------------------------- *

parameters
*   generator characteristics -- there data come from Geninmonth.mat
        Pmin(g)              minimum power generation of generator g
        Pmax(g)              maximum power generation of generator g
        FuelType(g)          fuel type (1-20) of generator g
        HeatRate(g)          heat rate of generator g        
        Avail(g,d)           renewable availability of generator g in hour d of week w
        GenArea(g)           generator by economic region (g \in R)
        a(g)                 production cost (variable) for generator g        
        Startupcost(g)       startup cost of generator g
        FuelCost(g)          fuel cost of production (variable) for generator g
        VOM(g)               variable operating and maintenance cost for generator g 
        Ramp(g)              ramp rate of generator g
        Uptime(g)            uptime requirement of generator g 
        Downtime(g)          downtime requirement of generator g
        GenLevel(gOff,d)     generation outage of generator g from WBM
        
*   demand parameters
        p1Demand(i,d)        demand at bus i in hour d of week w (demand in 2010)
*        DChange(GenLocation) percentage change in demand from REM
        DChange(zones) percentage change in demand from REM for UC model
        DChange_DC(econr) percentage change in demand from REM for UC model

;

* ------------------------------------------------------------------------------
*		Define Alias
* ------------------------------------------------------------------------------ 

*       alias for bus set
alias (i,j);

*       alias for zones set
alias(zones,zones2);

* ---------------------------------------------------------------------------- *
*       import sets and parameters from ParallelizingTrial.m:
* ---------------------------------------------------------------------------- *
$if not set gdxin $set gdxin MtoG
$GDXIN %gdxin%
$LOAD  i,g,ft,zones,econr,gzones,geconr,gi,gftype,ig,Gentozone
$LOAD  p1Demand,Avail
$LOAD  Pmin,Pmax,FuelType,GenArea,a,Startupcost,HeatRate,FuelCost,VOM,Ramp,Uptime,Downtime
$LOAD  gOff,GenLevel
$LOAD  DChange,Dchange_DC
$GDXIN

    display g, p1Demand, Pmin, Pmax;

* ------------------------------------------------------------------------------
*       load transmission line characteristic data from a text file
* ------------------------------------------------------------------------------

$ontext
loads
    TABLE dtl(line,i,j,cl) -- INFO?
$offtext

$include WECC_line_data2.txt

* ------------------------------------------------------------------------------
*		load transmission lines capacities and geographical crosswalks for line/lines
* ------------------------------------------------------------------------------ 

$ontext
defines
    LineCap(lines)                zone transmission line capacity
    LineSources(lines,zones)      zone transmission line starting zones crosswalk
    LineSink(lines,zones)         zone transmission line ending zones crosswalk
    Bustozone(i,zones)            bus to zone crosswalk (this can come from the Geninmonth data)
    linetozone(line,zones,zones2) zone transmission line connection
$offtext

$include trialimport_zline.gms

* ------------------------------------------------------------------------------
*       load bus to economic region crosswalk
* ------------------------------------------------------------------------------

$ontext
defines
    ieconr(i,econr) bus to economic region crosswalk
    BusZone(i)      bus i in economic region
    zeconr(zones,econr) zone zones to economic region crosswalk -- add this!!
$offtext

$include trialimport_econr.gms

* ------------------------------------------------------------------------------
*       load retired generator data from a gms file
* ------------------------------------------------------------------------------

$include trialimport_gretire.gms

* ------------------------------------------------------------------------------
*		define demand of zone zones for unit commitment program
* ------------------------------------------------------------------------------ 

parameters 

		LoadDemand(zones,d)                 demand of zone zones in hour d of week w
;

LoadDemand(zones,d)=sum(i$Bustozone(i,zones),p1Demand(i,d));

    display LoadDemand;


* ------------------------------------------------------------------------------
*		Re-define parameters
* ------------------------------------------------------------------------------  

*       for non-renewable generators, AvailFinal will be set to 1 by dividing Avail(g,d) by 10000
parameters
        AvailFinal(g,d) renewable availability of generator g in hour d of week w - removes 10000 
;

AvailFinal(g,d)=Avail(g,d)/10000;

    display AvailFinal;


* ------------------------------------------------------------------------------
*		define different generator dynamic sets
* ------------------------------------------------------------------------------ 

sets
        gNuke(g)   all nuclear generators
            /
                g98    Arizona	    Nuclear
                g99    Arizona	    Nuclear
                g100   Arizona	    Nuclear
                g1624  California	Nuclear
                g1625  California	Nuclear
                g1829  Washington	Nuclear
            /
            
        gg(g)      dynamic set for nuclear and coal retirement generators
        ggg(g)     dynamic set for california (GenArea==1) coal (FuelType==1) generators
        gggg(g)    "dynamic set for coal, gas CT, gas CCGT, steam turbine, biogas-other, petroleum and geothermal"
;

*       define dynamic set of california coal generators
ggg(g)$([(GenArea(g)=1) and (FuelType(g)=1)])=yes;

    display ggg;

* ------------------------------------------------------------------------------
*		define which shocks to turn on and off
* ------------------------------------------------------------------------------

*       baseline specification

*       UNCOMMENT to ENABLE RETIREMENTS
gg(g)$[gNuke(g) or gRetire(g)]=yes;

*       COMMENT to turn ON retirements
*gg(g)=no;


*       UNCOMMENT to turn OFF water shock
GenLevel(gOff,d) = 0;


    display gg;

* ------------------------------------------------------------------------------
*		define dynamic set of polluting generators 
* ------------------------------------------------------------------------------
   
*               add all generators into the polluting set
gggg(g) = yes;

*               remove renewables from set
$ontext
keep:
        Coal-PC == 1;
        Gas CT == 5;        Gas CCGT == 6
        Steam Turbine == 7; Biogas-Other == 14
        Petroleum == 16;    Geothermal == 17
        

remove:
        Wind-Onshore == 2; DG-BTM == 3; Hydro   == 4;
        Solar PV-Distributed Utility-Fixed Tilt == 8;
        Solar PV-Distributed Utility-Tracking   == 9;
        DC-Intertie == 10; DR == 11; Nuclear   == 12;
        PS-Hydro == 13;    SolarThermal == 15;
        ES-Generic  == 18; EE == 19;    MotorLoad == 20
$offtext

gggg(g)$([(FuelType(g)=2) or (FuelType(g)=3)  or (FuelType(g)=4) or (FuelType(g)=8)  or (FuelType(g)=9)  or (FuelType(g)=10) or (FuelType(g)=11)])=no;
gggg(g)$([(FuelType(g)=12)or (FuelType(g)=13) or (FuelType(g)=15)or (FuelType(g)=18) or (FuelType(g)=19) or (FuelType(g)=20) ])=no;

*               remove the retirements set from the polluting set
gggg(gg)=no;

    display gggg;
