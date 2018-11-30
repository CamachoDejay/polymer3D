function [pos] = nMaxDetection (data,nMax,minDist)

pos = [];
str_em = strel('disk', minDist);
while size(pos,1)< nMax
    
   max2Find = nMax - size(pos,1);
   maximum = max(max(data));
   [row,col] = find(data==maximum);
   
   if length(row) >1
       
       row = row(1);
       col = col(1);
       
   end
   
   pos = [pos; row col];
   BW = false(size(data));
   BW(row,col) = true;
   BW = imdilate(BW,str_em);
   data(BW==1) = 0;  
  
end