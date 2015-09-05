function [P_C] = pro_cpx(prot)
	% The function generate 2 order protein complex.
	% Example: input (y1,y2),output (1,y1,y2,y1*y2)
	% input(y1,y2,y3),output(1,y1,y2,y3,y1*y2,y1*y3,y2*y3)
    %%% Generate P_I up to second order

    [r,c] = size(prot); % number of proteins
    P_C = zeros((r*(r+1)/2)+1,c); % (c*(c+1)/2)+1 is the number of protein products  
    P_C(1,:) = 1; % no protein 
    
    %%% first order
    P_C(2:r+1,:) = prot; % single protein, no product
    
    %%% second order
    count = 1;
    for i = 1 : r-1
        for j = i+1 : r
            P_C(r+1+count,:) = prot(i,:).*prot(j,:);
            count = count+1;
        end
    end
end