classdef isochromat
    %This class defines a isochromat within the rotating FoR
    %   An individual isochromat with settable properties and a suite of
    %   methods that can simulate pulse-sequences
    
    properties (SetAccess=public, GetAccess=public)
        T2 = 1;
        T1 = 2;
        B0 = [0, 0, 1];
        x = [];
        y = [];
    end
    
    properties (SetAccess=private, GetAccess=public)
        Time = [0];
        Mx = [0];
        My = [0];
        Mz = [1];
        tRelax = [0];
    end
    
    methods
        function obj = InitialDir(obj, Axis)
            %Initialise the spin direction
            %   User defined direction of the spin, can only be set to
            %   along the positive axis directions
            if (nargin ~= 2)
                error('The number of arguments is: %d there should be only one argument either x, y or z.', nargin - 1)
            end
            if any((Axis~='x')&(Axis~='y')&(Axis~='z'))
                error('Needs be an axis of either "x","y" or "z".')
            end
            if (Axis=='x')
                M = [1, 0, 0];
            elseif (Axis=='y')
                M = [0, 1, 0];
            elseif (Axis=='z')
                M = [0, 0, 1];
            end
            obj.Mx = M(1);
            obj.My = M(2);
            obj.Mz = M(3);
        end
        
        function obj = Flip(obj, Axis, Angle)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if any((Axis~='x')&(Axis~='y')&(Axis~='z'))
                error('Needs be an axis of either "x","y" or "z".')
            end
            t = linspace(0,obj.FlipTime,obj.N);
            if (Axis=='x')
                M_x = obj.Mx(end)*ones(1,obj.N);
                M_y = obj.My(end).*cos(Angle/obj.FlipTime*t) + obj.Mz(end).*sin(Angle/obj.FlipTime*t);
                M_z = obj.Mz(end).*cos(Angle/obj.FlipTime*t) - obj.My(end).*sin(Angle/obj.FlipTime*t);
            elseif (Axis=='y')
                M_x = obj.Mx(end).*cos(Angle/obj.FlipTime*t) - obj.Mz(end).*sin(Angle/obj.FlipTime*t);
                M_y = obj.My(end).*ones(1,obj.N);
                M_z = obj.Mz(end).*cos(Angle/obj.FlipTime*t) + obj.Mx(end).*sin(Angle/obj.FlipTime*t);
            elseif (Axis=='z')
                M_x = obj.Mx(end).*cos(Angle/obj.FlipTime*t) + obj.My(end).*sin(Angle/obj.FlipTime*t);
                M_y = obj.My(end).*cos(Angle/obj.FlipTime*t) - obj.Mx(end).*sin(Angle/obj.FlipTime*t);
                M_z = obj.Mz(end)*ones(1,obj.N);
            end
            obj.Mx = [obj.Mx, M_x];
            obj.My = [obj.My, M_y];
            obj.Mz = [obj.Mz, M_z];
            obj.Time = [obj.Time, obj.FlipTime];
        end
        
        function obj = FreeRelaxZ(obj, omegaFoR)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if (obj.B0(3)~=1)
                error('B0 needs to be a unit vector along z-axis.')
            end
            omega_0 = omegaFoR - obj.omega0;
            t_Relax = linspace(sum(obj.tRelax),relaxTime+sum(obj.tRelax(end)),obj.N);
            
            M_z_begin = (obj.Mz(end)-(1-exp(-t_Relax(1)./obj.T1)))/exp(-t_Relax(1)/obj.T1);
            t = linspace(0,relaxTime,obj.N);
            M_x = (obj.Mx(end).*cos(omega_0*t) + obj.My(end).*sin(omega_0*t));%.*exp(-t_Relax/obj.T2);
            M_y = (obj.My(end).*cos(omega_0*t) - obj.Mx(end).*sin(omega_0*t));%.*exp(-t_Relax/obj.T2);
            M_z = M_z_begin*exp(-t_Relax/obj.T1) + (1-exp(-t_Relax./obj.T1));
            
            obj.Mx = [obj.Mx, M_x];
            obj.My = [obj.My, M_y];
            obj.Mz = [obj.Mz, M_z];
            obj.Time = [obj.Time, relaxTime];
            obj.tRelax = [obj.tRelax, relaxTime];
        end
        
        %         function M_x, M_y, M_z = method1(obj,[Mx,My,Mz],frequency, time)
        %             %METHOD1 Summary of this method goes here
        %             %   Detailed explanation goes here
        %             M_x = Mx*ones(1,obj.N);
        %             M_y = My.*cos(frequency*time) + obj.Mz(end).*sin(Angle/obj.FlipTime*t);
        %             M_z = Mz.*cos(Angle/obj.FlipTime*t) - obj.My(end).*sin(Angle/obj.FlipTime*t);
        %         end
        
        
        
    end
end

