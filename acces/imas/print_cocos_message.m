function  print_cocos_message(info,root)

noms = fieldnames(info);
for k = 1:length(noms)
    if isempty(info.(noms{k}))
        fprintf('Data %s is mentionned without any comment\n',noms{k})
    elseif iscell(info.(noms{k}))
        for l=1:length(info.(noms{k}))
            if isstruct(info.(noms{k}){l})
                    print_cocos_message(info.(noms{k}){l},sprintf('%s.%s{%d}',root,noms{k},l));
            else
                    fprintf('%s: %s\n',sprintf('%s.%s{%d}',root,noms{k},l),info.(noms{k}){l});               
            end
        end
    elseif isstruct(info.(noms{k}))
        print_cocos_message(info.(noms{k}),sprintf('%s.%s',root,noms{k}));
    else
        fprintf('%s: %s\n',sprintf('%s.%s',root,noms{k}),info.(noms{k}));
    end
end

