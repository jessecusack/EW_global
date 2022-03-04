function chbad = find_useful_moorings(meta_outfile,eps)

chbad = 0;


eval(['load ' meta_outfile]);


flag_C = zeros(size(Lon)); % checked
Num_moorings = 0; % number of moorings

Max_moorings = length(Lon);
% Max_mooring_inst = 1000;
% Mooring_instruments = NaN*ones(Max_moorings,Max_mooring_inst);
% Mooring_depths      = NaN*ones(Max_moorings,Max_mooring_inst);
% Mooring_lon         = NaN*ones(Max_moorings,1);
% Mooring_lat         = NaN*ones(Max_moorings,1);
% Mooring_overlap     = NaN*ones(Max_moorings,1);
% Mooring_ninst       = NaN*ones(Max_moorings,1);

[N,~,T,M,I] = get_metadata_list;
q = find(T==2);
for i = 1:length(q)
    if isequal(class(I{q(i)}),'double')
        Max_mooring_inst = I{q(i)};
        eval([N{q(i)}, ' = ', num2str(I{q(i)}), '*ones(', ...
            num2str(Max_moorings), ',', num2str(M(q(i))), ');']);
    else
        error('unexpected initialization class.')
    end
end

max_insts_on_mooring = 1;

Inst_mooring_num = NaN*Lon; % specifies which mooring each instrument belongs to
for j=1:length(Lon) % loop over instruments
    if flag_C(j)<1 % not yet checked
        

        nd0 = find(flag_C<1);
        % latitude difference:
        dy = abs(Lat(nd0)-Lat(j));
        % longitude difference:
        dx = abs(Lon(nd0)-Lon(j)); 
        nd1 = find(dx>180); 
        if length(nd1)>0; 
            dx(nd1) = 360-dx(nd1); 
        end;
        % differences in start and end time:
        dt1= abs(BeginTime(nd0)-BeginTime(j));
        dt2= abs(EndTime(nd0)-EndTime(j));
        % which instruments are at the same time and (horizonal) location?
        % Note: not checking depth this time; it is expected that multiple
        % instruments at a single mooring will have different depths.
        redund = find(dy < eps.xy & dx < eps.xy & dt1 < eps.t & dt2 < eps.t);
        % Hopefully, there is at least one instrument at this mooring.
        if length(redund)<1; disp(j)
            disp('Data not found')
            chbad = 1; 
            return
        end
        % Note: This error should never happen.
        flag_C(nd0(redund)) = 1; % mark these instruments as having already been grouped into moorings.
        Num_moorings = Num_moorings + 1; % keep track of how many moorings you have.
        Inst_mooring_num(nd0(redund)) = Num_moorings; % Tag all instruments in the mooring with the mooring number.
        % ch = NaN*ones(1,60); % Initial blank vector. Assumes that there are 60 or fewer instruments in all moorings.
        % ch(1:length(redund)) = nd0(redund); % Keep track of instrument numbers to associate with this mooring.
        % Mooring_instruments(Num_moorings,:) = ch;
        % ch(1:length(redund)) = Depth(nd0(redund)); % Keep track of depths of instruments in mooring.
        % Mooring_depths(Num_moorings,:) = ch;
        % ch = NaN*ones(1,Max_mooring_inst); % Initial blank vector. Assumes that there are 60 or fewer instruments in all moorings.
        % ch(1:length(redund)) = nd0(redund); % Keep track of instrument numbers to associate with this mooring.
        Mooring_instruments(Num_moorings,1:length(redund)) = nd0(redund)';
        % ch(1:length(redund)) = Depth(nd0(redund)); % Keep track of depths of instruments in mooring.
        Mooring_depths(Num_moorings,1:length(redund)) = Depth(nd0(redund));
        Mooring_lat(Num_moorings) = mean(Lat(nd0(redund)));
        Mooring_lon(Num_moorings) = mean(Lon(nd0(redund)));
        Mooring_overlap(Num_moorings) = min(EndTime(nd0(redund)))-max(BeginTime(nd0(redund)));
        Mooring_ninst(Num_moorings) = length(redund);
        
        if length(redund) > max_insts_on_mooring
            max_insts_on_mooring = length(redund);
        end
        
        
    end
end % end loop over instruments

if max(Mooring_ninst) > max_insts_on_mooring
    error('Maximum number of instruments on a mooring under estimated.')
end

% Mooring_instruments = ...
%     Mooring_instruments(1:Num_moorings,1:max(Mooring_ninst));
% Mooring_depths      = Mooring_depths(1:Num_moorings,1:max(Mooring_ninst));
% Mooring_lat = Mooring_lat(1:Num_moorings,1);
% Mooring_lon = Mooring_lon(1:Num_moorings,1);
% Mooring_ninst = Mooring_ninst(1:Num_moorings,1);

[N,~,T,M] = get_metadata_list;
q = find(T==2);
for i = 1:length(q)
    if M(q(i)) > 1
        eval([N{q(i)}, '=', N{q(i)}, '(1:', num2str(Num_moorings), ',', ...
            '1:', num2str(max_insts_on_mooring), ');']);
    else
        eval([N{q(i)}, '=', N{q(i)}, '(1:', num2str(Num_moorings), ...
            ',1);']);
    end
end



% add mooring information to the user's output metadata file:
% eval(['save ' meta_outfile ' Mooring_lat Mooring_lon ',...
%     'Num_moorings Inst_mooring_num Mooring_instruments',...
%     ' Mooring_depths Mooring_ninst -append']);
[N,~,T] = get_metadata_list;
q = find(T==2 | T==3);
for i = 1:length(q)
   eval(['save  ', meta_outfile, ' -append ', N{q(i)}]);
end

