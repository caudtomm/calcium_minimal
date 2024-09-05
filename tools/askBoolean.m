function outbit = askBoolean(prompt)

while true
    x = input(prompt, 's');
    if strcmp(x,'y') || strcmp(x,'Y') || strcmp(x,'1') || strcmp(x,'yes')
        outbit = true; break
    elseif strcmp(x,'n') || strcmp(x,'N') || strcmp(x,'0') || strcmp(x,'no')
        outbit = false; break
    else
        disp('Input invalid!')
    end
end