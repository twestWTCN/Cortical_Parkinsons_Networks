function [amp phi dphi_12 dphi_12_dt] = comp_instant_angle_phase_nobp(x)

        
        amp(:,1) = abs(hilbert(x(1,:)));
        amp(:,2) = abs(hilbert(x(2,:)));
         amp(:,2) = abs(hilbert(x(2,:)));
        phi(:,1) = angle(hilbert(x(1,:)));
        phi(:,2) = angle(hilbert(x(2,:)));
        phi(:,3) = angle(hilbert(x(2,:)));
        dphi_12 = unwrap(phi(:,1)-phi(:,2));
        % optional amp weighting
%         ampw = unwrap((amp(:,1).*amp(:,2))./(max(amp(:,1))*max(amp(:,2))));  %%
%         dphi_12 = angle(ampw.*exp(1i.*(phi(:,1)-phi(:,2))));%%
%         dphi_12 = (dphi_12-(1/sqrt(length(phi(:,1)))))./(1-(1/sqrt(length(phi(:,1))))); %%
        dphi_12_dt = diff(dphi_12);
        
