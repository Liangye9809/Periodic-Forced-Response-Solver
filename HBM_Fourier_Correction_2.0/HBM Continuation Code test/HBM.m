classdef HBM
    properties
        H
        N
        Nc
        Nx
        Na
        xi
        E
        EH
        Idx
        H_F_ext
        fftfa
        fftfx
    end

    methods
        function obj = HBM(HBMstruct)
            obj.H = HBMstruct.H;
            obj.N = HBMstruct.N;
            obj.Nc = HBMstruct.Nx; % dof of contact part
            obj.Nx = HBMstruct.Nx; % dof of contact part
            obj.Na = HBMstruct.Na;
            obj.xi = HBMstruct.xi;
            obj.H_F_ext = HBMstruct.H_F_ext;
            obj.Idx = HBM.Reorder(HBMstruct.H, 3 * HBMstruct.Nx, HBMstruct.Na);  % Call static method
            [obj.E, obj.EH] = obj.fft_matrices(HBMstruct.N, HBMstruct.H);
            [obj.fftfa, obj.fftfx] = HBM.fftf(HBMstruct.H, HBMstruct.H_F_ext, HBMstruct.CB_F);
        end
    end

    methods (Static, Access = private)
        function idx = Reorder(H, Nx, Na)
            n = 2 * H + 1;

            Ia = (1:n * Na)';
            Idxa = Ia(1:Na:end);
            for i = 2:Na
                Idxa = [Idxa; Ia(i:Na:end)];
            end

            Iax = (1:n * Nx)';
            Idxx = Iax(1:Nx:end);
            for i = 2:Nx
                Idxx = [Idxx; Iax(i:Nx:end)];
            end

            idx = [Idxa; Idxx + n * Na];
        end

        function [fftfa, fftfx] = fftf(H, H_F_ext, CB_F)
            H_F_ext = [H_F_ext, zeros(1, 2*H+1 - size(H_F_ext,2))];
            Fa = CB_F.Fa * H_F_ext;
            Fa = Fa';
            fftfa = Fa(:); 
            Fx = CB_F.Fx * H_F_ext;
            Fx = Fx';
            fftfx = Fx(:);
        end

    end

    methods (Static)
        function [E, EH] = fft_matrices(N, H)
            E = zeros(N, 2 * H + 1);
            E(:, 1) = 1 / 2;
        
            t = (0:(N - 1))';
            
            for h = 1:H
                E(:, 2 * h) = cos(2*pi/N * t * h);
                E(:, 2 * h + 1) = sin(2*pi/N * t * h);
            end
            
            EH = E';
            EH(1,:) = 1;
            
            EH = EH * (2/N);
        end
    end
end
