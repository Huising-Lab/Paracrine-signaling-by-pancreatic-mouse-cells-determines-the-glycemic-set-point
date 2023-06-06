% Created by: Mohammad Pourhosseinzadeh
% Huising Lab

% This function is used to analyze single cell Ca2+ imaging data. It finds the first
% instance at which the fluorescence intensity of a single cell passes half maximal intensity,
% typically defined by a brief pulse of KCl.

%------------------------------------------------------------------------------%
% INPUT: 
% file= Name of the csv file minus the file extension
%       -It is expected that the first column is the a column vector containing the time of each measurement in seconds
%       -Every column after the first is expected to be raw Ca2+ data
%       -Expected to have no headers in the data file
% t1 and t2= Upper and lower bounds of window for data analysis
%------------------------------------------------------------------------------%

function [hlf_mx]=hlf_max(file,t1,t2)
  data_import=dlmread(strcat(file,'.csv'));
  time=data_import(:,1);

  t1=t1.*60;
  t2=t2.*60;
  
  % Find the index of t1 and t2 in the vector 'time'
  [val1,trow1]=min(abs(time-t1));
  [val2,trow2]=min(abs(time-t2));
  
  data=data_import(trow1:trow2,2:end);
  time=time(trow1:trow2,1);
   
  hlf_mx=[];
  i=1;
  while i<=columns(data);
    ndata=(data(:,i)-min(data(:,i)))./(max(data(:,i)-min(data(:,i)))); % Normalize the data
    temp=ndata-0.5; % subtract the threshold 0.5
    temp_prod=temp.*[-1;temp(1:rows(ndata)-1,1)]; %multiply the data by itself shifted down by one index 
    [indx]=find(temp_prod<=0); % find the index of all negative points in temp_prod as those are all the points that have passed the threshold 0.5
    hlf_mx(i,1)=time(indx(1,1),1)-t1; % Convert the index of these points into time (s)
    i++; 
  endwhile
  csvwrite(strcat(file,'_half_max','.csv'),[hlf_mx])
endfunction