classdef CoulombFrictionParas
    
    properties

        kn
        xn0
        mu
        kt
        w

    end

    methods
        function obj = CoulombFrictionParas(Coulombstruct)

            Nx = Coulombstruct.Nx;
            kn = Coulombstruct.kn;
            xn0 = Coulombstruct.xn0;
            mu = Coulombstruct.mu;
            kt = Coulombstruct.kt;
            % w = [0;0];
            w = Coulombstruct.w;

            if size(kn, 2) == 1

                obj.kn = kn * ones(1, Nx);
                obj.xn0 = xn0 * ones(1, Nx);
                obj.mu = mu * ones(1, Nx);
                obj.kt = kt * ones(1, Nx);
                obj.w = w * ones(1, Nx);

            elseif size(kn, 2) == Nx

                obj.kn = kn;
                obj.xn0 = xn0;
                obj.mu = mu;
                obj.kt = kt;
                obj.w = w;
                
            else
                error(['Number of coulomb friction parameters is not match with the contact point' ...
                    'when evey contact point has different of that']);
            end
            
        end
    end
end

