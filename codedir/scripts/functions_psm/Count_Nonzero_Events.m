function [event_data] = Count_Nonzero_Events(demanddata, nsedata, timeslice32)

event_data = cell(32,8);

DATALEN = length(demanddata);

EVENT_FLAG = 0;
counthrs = 0;
demandsum = 0;
nsesum = 0;

for t = 1:DATALEN
   if  (nsedata(t))
       % non-zero value, check if already counting
       if (EVENT_FLAG)
           % continue counting
           counthrs = counthrs + 1;
           demandsum = demandsum + demanddata(t);
           nsesum = nsesum + nsedata(t);
       else
           % start counting new event
           counthrs = 1;
           demandsum = demanddata(t);
           nsesum = nsedata(t);
           EVENT_FLAG = 1;
       end
   else
       % zero NSE value, check if just finished counting
       if (EVENT_FLAG)
           % just finished an event, so now store and reset
           avg_nse = nsesum / demandsum;
           
           % FOR NOW allocate event to timeslice of final hour
           event_time = timeslice32(t-1);
           
           % how long was the event?
           if (counthrs > 8)
               event_length = 8;
           else
               event_length = counthrs;
           end
           
           event_data{event_time, event_length} = [event_data{event_time, event_length} avg_nse];
           
           % reset for next event
           counthrs = 0;
           demandsum = 0;
           nsesum = 0;
           EVENT_FLAG = 0;
       else
           % do nothing, continuing run of zeros
       end
   end
   
           
           
    
end




end

