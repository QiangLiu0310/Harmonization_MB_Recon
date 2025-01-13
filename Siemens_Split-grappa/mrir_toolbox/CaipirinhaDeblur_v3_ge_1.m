function [Kcorrected] = CaipirinhaDeblur_v3_ge_1(K)

PhaseShiftPerMM = (2*pi/3/90); 

SlicePos=[44.5:-1:0];
SlicePos  = SlicePos(end:-1:1); 
PhaseShiftBtwSimulSlices= 2*pi/3 ;

Kcorrected = (zeros(size(K)));

% if PhaseShiftBtwSimulSlices ~= 0
    for SlcCount = 1: size(K,3), % ( length(prot.sSliceArray.asSlice) > 1 ),
      PhaseShift = PhaseShiftPerMM*SlicePos(SlcCount);
        if abs(PhaseShiftBtwSimulSlices) == pi % FOV/2 shift
            Kcorrected(:,1:2:end,SlcCount,:) =  K(:,1:2:end,SlcCount,:);
            Kcorrected(:,2:2:end,SlcCount,:) =  K(:,2:2:end,SlcCount,:)*exp(-i*PhaseShift);
         elseif abs(PhaseShiftBtwSimulSlices) == 2*pi/3 % FOV/3 shift
            Kcorrected(:,1:3:end,SlcCount,:) =  K(:,1:3:end,SlcCount,:);
            Kcorrected(:,2:3:end,SlcCount,:) =  K(:,2:3:end,SlcCount,:)*exp(-i*PhaseShift);
            Kcorrected(:,3:3:end,SlcCount,:) =  K(:,3:3:end,SlcCount,:)*exp(-i*2*PhaseShift);
        elseif abs(PhaseShiftBtwSimulSlices) == pi/2 % FOV/4 shift
            Kcorrected(:,1:4:end,SlcCount,:) =  K(:,1:4:end,SlcCount,:);
            Kcorrected(:,2:4:end,SlcCount,:) =  K(:,2:4:end,SlcCount,:)*exp(-i*PhaseShift);
            Kcorrected(:,3:4:end,SlcCount,:) =  K(:,3:4:end,SlcCount,:)*exp(-i*2*PhaseShift);
            Kcorrected(:,4:4:end,SlcCount,:) =  K(:,4:4:end,SlcCount,:)*exp(-i*3*PhaseShift);
        end
    end



