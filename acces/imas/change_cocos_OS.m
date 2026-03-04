function [ids_out,cocoscoeff,cocoscoeffsym] = change_cocos_OS(ids_in,ids_type,cocos_struct)

% for DDv4 swap
switch ids_type
    case 'plasma_profiles'
        ids_type = 'core_profiles';
    case 'plasma_sources'
        ids_type = 'core_sources';
end
%
% empty result is not call
ids_cocos = [];
stat_struct  = [];
end_verif = false;

if cocos_struct.COCOS_Sauter
    % test for supported IDS
    %        
    switch ids_type
        case {'equilibrium', 'magnetics', 'pf_active','wall', 'tf','core_profiles','core_sources','ec_launchers','nbi','pf_passive','summary'}
            fprintf('COCOS mapping for %s\n',ids_type)
            [ids_out,cocoscoeff]=ids_generic_cocos_nodes_transformation_symbolic(ids_in,ids_type,cocos_struct.cocos_in,cocos_struct.cocos_out,cocos_struct.ipsign_out, ...
                cocos_struct.b0sign_out,cocos_struct.ipsign_in,cocos_struct.b0sign_in,cocos_struct.error_bar,cocos_struct.verbose);
            % some field are not handle by the tool for equilibrium
            switch ids_type
                case 'equilibrium'
                    for lm = 1:length(ids_out.time_slice)
                        sbt = sign(ids_out.vacuum_toroidal_field.b0(lm));
                        for ln=1:length(ids_out.time_slice{lm}.profiles_2d)
                            if isfield(ids_out.time_slice{lm}.profiles_2d{ln},'b_tor')
                                ids_out.time_slice{lm}.profiles_2d{ln}.b_tor = sbt .* sign(mean(ids_out.time_slice{lm}.profiles_2d{ln}.b_tor(:))) .* ids_out.time_slice{lm}.profiles_2d{ln}.b_tor;
                            end
                            if isfield(ids_out.time_slice{lm}.profiles_2d{ln},'b_phi')
                                ids_out.time_slice{lm}.profiles_2d{ln}.b_phi = sbt .* sign(mean(ids_out.time_slice{lm}.profiles_2d{ln}.b_phi(:))) .* ids_out.time_slice{lm}.profiles_2d{ln}.b_phi;
                            end
                            if isfield(ids_out.time_slice{lm}.profiles_2d{ln},'b_field_tor')
                                ids_out.time_slice{lm}.profiles_2d{ln}.b_tor = sbt .* sign(mean(ids_out.time_slice{lm}.profiles_2d{ln}.b_field_tor(:))) .* ids_out.time_slice{lm}.profiles_2d{ln}.b_field_tor;
                            end
                            if isfield(ids_out.time_slice{lm}.profiles_2d{ln},'b_field_phi')
                                ids_out.time_slice{lm}.profiles_2d{ln}.b_field_phi = sbt .* sign(mean(ids_out.time_slice{lm}.profiles_2d{ln}.b_field_tor(:))) .* ids_out.time_slice{lm}.profiles_2d{ln}.b_field_phi;
                            end
                        end
                    end
            end
            
            switch  cocos_struct.check
                case 'on'
                    % multi dictionnary version compatibility
                    try
                        % backward compatibility
                        model = ids_gen(ids_type);
                    catch
                        model = [];
                    end
                    if ~isempty(model)
                        model_fields = fieldnames(model);
                        ids_out_clear = ids_out;
                        data_fields  = fieldnames(ids_out_clear);
                        for lk =1:length(data_fields)
                            if isempty(strmatch(data_fields{lk},model_fields,'exact'))
                                fprintf('Non existing field %s in IMAS model\n',data_fields{lk});
                                ids_out_clear = rmfield(ids_out_clear,data_fields{lk});
                            end
                        end
                    end
                    % start test
                    switch ids_type
                        case {'core_profiles','equilibrium'}
                            disp('****************************************************')
                            disp('METIS is natively in COCOS = 7 with ipsign = 1 and b0sign = 1');
                            disp('COCOS transform parameters:');
                            cocos_struct
                            disp('COCOS transformation coeff are:')
                            cocoscoeff
                    end
                    switch ids_type
                        case 'core_profiles'
                            disp('---------------------------------------')
                            disp('Verifiing COCOS in IDS core_profiles');
                            [ids_cocos,cocos_struct_out,stat_struct]=ids_generic_cocos_check(ids_out_clear,ids_type,cocos_struct.cocos_out,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
                            end_verif = true;
                        case 'equilibrium'
                            disp('---------------------------------------')
                            disp('Verifiing COCOS in IDS iequilibrium');
                            % keep only inverse grid others are not tested
                            for l=1:length(ids_out_clear.time_slice)
                                ids_out_clear.time_slice{l}.profiles_2d = ids_out_clear.time_slice{l}.profiles_2d(2);
                            end
                            [ids_cocos,cocos_struct_out,stat_struct]=ids_generic_cocos_check(ids_out_clear,ids_type,cocos_struct.cocos_out,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
                            end_verif = true;
                   end
             end
            
            if nargout > 2
                [~,cocoscoeffsym]=cocos_transform_coefficients(cocos_struct.cocos_in, cocos_struct.cocos_out,[],[], ...
                    cocos_struct.ipsign_in,cocos_struct.b0sign_in,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
            end
        otherwise
            ids_out = ids_in;
            if nargout > 1
                [cocoscoeff,cocoscoeffsym]=cocos_transform_coefficients(cocos_struct.cocos_in, cocos_struct.cocos_out,[],[], ...
                    cocos_struct.ipsign_in,cocos_struct.b0sign_in,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
            end
    end
    
else
    ids_out = ids_in;
    if (nargout > 1) || strcmp(cocos_struct.check,'on')
        [cocoscoeff,cocoscoeffsym]=cocos_transform_coefficients(cocos_struct.cocos_in, cocos_struct.cocos_out,[],[], ...
                    cocos_struct.ipsign_in,cocos_struct.b0sign_in,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
    end
    
    switch  cocos_struct.check
        case 'on'
            % multi dictionnary version compatibility
            try
                % backward compatibility
                model = ids_gen(ids_type);
            catch
                model = [];
            end
            if ~isempty(model)
                model_fields = fieldnames(model);
                ids_out_clear = ids_out;
                data_fields  = fieldnames(ids_out_clear);
                for lk =1:length(data_fields)
                    if isempty(strmatch(data_fields{lk},model_fields,'exact'))
                        fprintf('Non existing field %s in IMAS model\n',data_fields{lk});
                        ids_out_clear = rmfield(ids_out_clear,data_fields{lk});
                    end
                end
            end
            % start test
            switch ids_type
                case {'core_profiles','equilibrium'}
                    disp('****************************************************')
                    disp('METIS is natively in COCOS = 7 with ipsign = 1 and b0sign = 1');
                    disp('COCOS transform parameters:');
                    cocos_struct
                    disp('COCOS transformation coeff are:')
                    cocoscoeff
            end
            switch ids_type
                case 'core_profiles'
                    disp('---------------------------------------')
                    disp('Verifiing COCOS in IDS core_profiles');
                    [ids_cocos,cocos_struct_out,stat_struct]=ids_generic_cocos_check(ids_out_clear,ids_type,cocos_struct.cocos_out,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
                    end_verif = true;
                case 'equilibrium'
                    disp('---------------------------------------')
                    disp('Verifiing COCOS in IDS iequilibrium');
                    [ids_cocos,cocos_struct_out,stat_struct]=ids_generic_cocos_check(ids_out_clear,ids_type,cocos_struct.cocos_out,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
                    end_verif = true;
           end
    end
end

if ~isempty(stat_struct)
    noms = fieldnames(stat_struct);
    for k = 1:length(noms)
        switch noms{k}
            case 'global_status'
                fprintf('Number of error = %d\n',stat_struct.global_status);
            otherwise
                print_cocos_message(stat_struct.(noms{k}),noms{k});
        end
    end
%     if stat_struct.global_status > 0
%         keyboard
%     end
end
if end_verif
    disp('End IDS verification')
    disp('****************************************************')
    disp(' ');
end


% switch ids_type
%     case 'equilibrium'
%         filename = sprintf('equilibrium_cin_%d_cout_%d_ipsin_%d_bspin_%d_ipso_%d_bso_%d.mat',cocos_struct.cocos_in, cocos_struct.cocos_out, ...
%                    cocos_struct.ipsign_in,cocos_struct.b0sign_in,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
%         save(filename,'ids_in','ids_out');
% end






















