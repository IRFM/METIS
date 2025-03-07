function [ids_out,cocoscoeff,cocoscoeffsym] = change_cocos_OS(ids_in,ids_type,cocos_struct)


if cocos_struct.COCOS_Sauter
    % test for supported IDS
    switch ids_type
        case {'equilibrium', 'magnetics', 'pf_active','wall', 'tf','core_profiles','core_sources','ec_launchers','nbi','pf_passive','summary'}
            fprintf('COCOS mapping for %s\n',ids_type)
            [ids_out,cocoscoeff]=ids_generic_cocos_nodes_transformation_symbolic(ids_in,ids_type,cocos_struct.cocos_in,cocos_struct.cocos_out,cocos_struct.ipsign_out, ...
                cocos_struct.b0sign_out,cocos_struct.ipsign_in,cocos_struct.b0sign_in,cocos_struct.error_bar,cocos_struct.verbose);
            switch  cocos_struct.check
                case 'on'
                    % multi dictionnary version compatibility
                    model = ids_gen(ids_type);
                    model_fields = fieldnames(model);
                    ids_out_clear = ids_out;
                    data_fields  = fieldnames(ids_out_clear);
                    for lk =1:length(data_fields)
                        if isempty(strmatch(data_fields{lk},model_fields,'exact'))
                            fprintf('Non existing field %s in IMAS model\n',data_fields{lk});
                            ids_out_clear = rmfield(ids_out_clear,data_fields{lk});
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
                            ids_generic_cocos_check(ids_out_clear,ids_type,cocos_struct.cocos_out,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
                        case 'equilibrium'
                            disp('---------------------------------------')
                            disp('Verifiing COCOS in IDS iequilibrium');
                            ids_generic_cocos_check(ids_out_clear,ids_type,cocos_struct.cocos_out,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
                   end
                    switch ids_type
                        case {'core_profiles','equilibrium'}
                            disp('End IDS verification')
                            disp('****************************************************')
                            disp(' ');
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
            model = ids_gen(ids_type);
            model_fields = fieldnames(model);
            ids_out_clear = ids_out;
            data_fields  = fieldnames(ids_out_clear);
            for lk =1:length(data_fields)
                if isempty(strmatch(data_fields{lk},model_fields,'exact'))
                    fprintf('Non existing field %s in IMAS model\n',data_fields{lk});
                    ids_out_clear = rmfield(ids_out_clear,data_fields{lk});
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
                    ids_generic_cocos_check(ids_out_clear,ids_type,cocos_struct.cocos_out,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
                case 'equilibrium'
                    disp('---------------------------------------')
                    disp('Verifiing COCOS in IDS iequilibrium');
                    ids_generic_cocos_check(ids_out_clear,ids_type,cocos_struct.cocos_out,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
            end
            switch ids_type
                case {'core_profiles','equilibrium'}
                    disp('End IDS verification')
                    disp('****************************************************')
                    disp(' ');
            end
    end

end

% switch ids_type
%     case 'equilibrium'
%         filename = sprintf('equilibrium_cin_%d_cout_%d_ipsin_%d_bspin_%d_ipso_%d_bso_%d.mat',cocos_struct.cocos_in, cocos_struct.cocos_out, ...
%                    cocos_struct.ipsign_in,cocos_struct.b0sign_in,cocos_struct.ipsign_out,cocos_struct.b0sign_out);
%         save(filename,'ids_in','ids_out');
% end