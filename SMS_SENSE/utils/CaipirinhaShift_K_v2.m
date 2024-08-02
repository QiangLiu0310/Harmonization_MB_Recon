function [Kcorrected] = CaipirinhaShift_K_v2(K, CurrentSliceGroup,PhaseShiftBase)

%% ________________________________________________________________________
%% Caipirinha for Slice GRAPPA with virtual coils
%% ________________________________________________________________________
%% 
%% Haifeng Wang 
%% Sept. 13 2016
%%
%% now work for all acquisition orientations
%% assume that the slices are already sorted (so can use para in original prot to calc stuff) 
%%
%% SliceSep in mm
%% ________________________________________________________________________

PhaseShift = rem(PhaseShiftBase*(CurrentSliceGroup), 2*pi);

Kcorrected = (zeros(size(K)));

if PhaseShift ~= 0 
    if abs(PhaseShiftBase) == pi % FOV/2 shift
        Kcorrected(:,1:2:end,:,:,:,:,:,:,:,:) =  K(:,1:2:end,:,:,:,:,:,:,:,:);
        Kcorrected(:,2:2:end,:,:,:,:,:,:,:,:) =  K(:,2:2:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
        %Kcorrected = Kcorrected*exp(-i*PhaseShift/2);
    elseif abs(PhaseShiftBase) == 2*pi/3 % FOV/3 shift
        Kcorrected(:,1:3:end,:,:,:,:,:,:,:,:) =  K(:,1:3:end,:,:,:,:,:,:,:,:);
        Kcorrected(:,2:3:end,:,:,:,:,:,:,:,:) =  K(:,2:3:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
        Kcorrected(:,3:3:end,:,:,:,:,:,:,:,:) =  K(:,3:3:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
        Kcorrected = Kcorrected*exp(-i*PhaseShift);
    elseif abs(PhaseShiftBase) == pi/2 % FOV/4 shift
        Kcorrected(:,1:4:end,:,:,:,:,:,:,:,:) =  K(:,1:4:end,:,:,:,:,:,:,:,:);
        Kcorrected(:,2:4:end,:,:,:,:,:,:,:,:) =  K(:,2:4:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
        Kcorrected(:,3:4:end,:,:,:,:,:,:,:,:) =  K(:,3:4:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
        Kcorrected(:,4:4:end,:,:,:,:,:,:,:,:) =  K(:,4:4:end,:,:,:,:,:,:,:,:)*exp(i*3*PhaseShift);
        %Kcorrected = Kcorrected*exp(-i*PhaseShift*(3/2));
    elseif abs(PhaseShiftBase) == pi/3 % FOV/6 shift
        Kcorrected(:,1:6:end,:,:,:,:,:,:,:,:) =  K(:,1:6:end,:,:,:,:,:,:,:,:);
        Kcorrected(:,2:6:end,:,:,:,:,:,:,:,:) =  K(:,2:6:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
        Kcorrected(:,3:6:end,:,:,:,:,:,:,:,:) =  K(:,3:6:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
        Kcorrected(:,4:6:end,:,:,:,:,:,:,:,:) =  K(:,4:6:end,:,:,:,:,:,:,:,:)*exp(i*3*PhaseShift); 
        Kcorrected(:,5:6:end,:,:,:,:,:,:,:,:) =  K(:,5:6:end,:,:,:,:,:,:,:,:)*exp(i*4*PhaseShift);
        Kcorrected(:,6:6:end,:,:,:,:,:,:,:,:) =  K(:,6:6:end,:,:,:,:,:,:,:,:)*exp(i*5*PhaseShift); 
        %Kcorrected = Kcorrected*exp(-i*PhaseShift*(5/2));
        Kcorrected = Kcorrected*exp(-i*PhaseShift);
    elseif abs(PhaseShiftBase) == pi/4 % FOV/8 shift    
        Kcorrected(:,1:8:end,:,:,:,:,:,:,:,:) =  K(:,1:8:end,:,:,:,:,:,:,:,:);
        Kcorrected(:,2:8:end,:,:,:,:,:,:,:,:) =  K(:,2:8:end,:,:,:,:,:,:,:,:)*exp(i*PhaseShift);
        Kcorrected(:,3:8:end,:,:,:,:,:,:,:,:) =  K(:,3:8:end,:,:,:,:,:,:,:,:)*exp(i*2*PhaseShift);
        Kcorrected(:,4:8:end,:,:,:,:,:,:,:,:) =  K(:,4:8:end,:,:,:,:,:,:,:,:)*exp(i*3*PhaseShift); 
        Kcorrected(:,5:8:end,:,:,:,:,:,:,:,:) =  K(:,5:8:end,:,:,:,:,:,:,:,:)*exp(i*4*PhaseShift);
        Kcorrected(:,6:8:end,:,:,:,:,:,:,:,:) =  K(:,6:8:end,:,:,:,:,:,:,:,:)*exp(i*5*PhaseShift); 
        Kcorrected(:,7:8:end,:,:,:,:,:,:,:,:) =  K(:,7:8:end,:,:,:,:,:,:,:,:)*exp(i*6*PhaseShift);
        Kcorrected(:,8:8:end,:,:,:,:,:,:,:,:) =  K(:,8:8:end,:,:,:,:,:,:,:,:)*exp(i*7*PhaseShift); 
        %Kcorrected = Kcorrected*exp(-i*PhaseShift*(7/2));
    else
        keyboard
    end       
else
    Kcorrected = K;
end
