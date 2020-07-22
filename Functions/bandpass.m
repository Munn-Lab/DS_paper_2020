
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[dataFiltered] = bandpass( data, lowFreq, highFreq, Fs )
%
% Filter data with FFT based bandpass. Bandpass frequencies are defined by
% 'lowFreq' and 'highFreq' (in Hz), and assuming the sampling frequency defined
% by Fs (in kHz).  Note from LLC on Jan 9 2006- I changed this to read the
% sampling frequency in Hz.
%
% data: each column is assumed to be a continuous trace.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dataFiltered] = bandpass( data, lowFreq, highFreq, Fs )


%--- reshape data to be vertical ---

if( size(data,1) < size(data,2) )
	data = data';
end 


%--- Filter parmeters ---

%Fs = Fs * 1000;		% convert to Hz
lengthTrace = size( data, 1 );
nCutLow  = (lowFreq * lengthTrace) / Fs + 1;
nCutHigh = (highFreq * lengthTrace) / Fs + 1;


%===== Filter the data =====

F = fft( data ); 

%fig(1)
%plot(abs(F))

cx0 = complex(0);

halfLength = floor( lengthTrace / 2 );


%--- low cut ---

if (1 < nCutLow)					% the zero component ...
   F( 1, : ) = cx0;          	    
end

for k = 2:nCutLow-1		
    F( k, : ) = cx0;    	      	    % Set this component to zero...
    F( lengthTrace - k + 2, : ) = cx0;	% ...for all channels
end


%--- high cut ---

for k = nCutHigh+1:halfLength+1			% the + 1 takes care of the even/odd cases
   F( round(k), : ) = cx0;                	    % Set this component to zero...
   F( lengthTrace - round(k) + 2, : ) = cx0;		% ...for all channels
end


%hold on, plot(abs(F),'r'), hold off


dataFiltered = real( ifft(F) );


