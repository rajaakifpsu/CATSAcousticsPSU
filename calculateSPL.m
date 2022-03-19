%Purpose: to calculate the SPL in a given frequency interval

%Inputs: freqs is vector of frequencies usually pwelch__ in other codes.
%       pxx1mic is pxx but only for 1 mic (ie only one column of the full
%       matrix variable elsewhere in code)
%       P_ref is reference pressure
%       interval is 2-member vector for left side and right limits for the
%       freq vector (interval is 2 indexes)

%Outputs: SPLreturn is a single value, the SPL of pxx1mic in the interval.



function SPLreturn=calculateSPL(freqs,pxx1mic,P_ref,interval)

%Doing reiman sum to aproximate the SPL of the tones corresponding to
%frequencies of BPF and its harmonics
    rsum=0;
    pre = 0;
    dfpwelch=freqs(2)-freqs(1);
    fprintf('Interval: %i to %i\n',interval(1),interval(2))
    for j=interval(1):(interval(2)-1)%loops through frequencies corresponding to the range desired (-1 because we do j+1 in loop)
        rsum=rsum+(.5 * (freqs(j+1)-freqs(j)) * (pxx1mic(j)+pxx1mic(j+1) ) );
        
        pre = pre + (dfpwelch*(pxx1mic(j) + pxx1mic(j+1)));
        rsum2 = pre;
    end
    
    %similar to OASPL calculation oaspl calculate 
    
    SPLreturn=10*log10(rsum/(P_ref^2)); 
    SPLreturn2=10*log10(rsum2/(P_ref^2)); 
 
    %hello
end