function [corrF,corrF2] = calcCorr(A,B,lag)
assert(length(A) == length(B), 'Vector should be the same length');
corrF = zeros(1,length(A)-lag);

idx = 1:lag:length(A);
for i=1:length(A)-lag
    
   corrF(i) = corr(A(i:i+lag),B(i:i+lag)); 
   
   if i< length(idx)
   startIdx = idx(i);
  
   if idx(i) == idx(end)
       endIdx = length(A);
   else
       endIdx = idx(i+1);
   end
       
   corrF2(startIdx:endIdx) = corr(A(startIdx:endIdx),B(startIdx:endIdx)); 
   else
   end
%     stp=i;cnt=1;tmpCorr=[];tmpCorr2=[];
% 
%     while cnt<=stp & cnt+stp<=length(A)
%         idx=cnt:stp:length(A);
%         for j = 1:length(idx)-1
%             tmpCorr2 = [tmpCorr2 corr(A(idx(j):idx(j+1)),B(idx(j):idx(j+1)))];
%             
%         end
%         
%         
%         tmpCorr =[tmpCorr corr(A(idx),B(idx))];
%         cnt=cnt+1;
%         
%     end %while
%     corrF(i) = mean(tmpCorr);
%     corrF2(i) = mean(tmpCorr2);
end



end