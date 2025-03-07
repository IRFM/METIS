liste          = {};
liste{end+1}   = 'ppf/@shot/BOLO/TOPI';    % puissance rayonnee totale improved
liste{end+1}   = 'ppf/@shot/BOLO/TOPO';    % puissance rayonnee totale old
liste{end+1}   = 'ppf/@shot/BOLO/TOBU';    % puissance rayonnee coeur  upper
liste{end+1}   = 'ppf/@shot/BOLO/TOBH';    % puissance rayonnee coeur  horizontal
liste{end+1}   = 'ppf/@shot/BOLO/TXPN';    % puissance rayonnee xpoint and divertor   new
liste{end+1}   = 'ppf/@shot/BOLO/TOXP';    % puissance rayonnee xpoint and divertor old
bolo         = cgcgetjet(post.z0dinput.shot,liste,'','');

