function arrayout = nanzscore(arrayin,dim)
    arguments
        arrayin double
        dim double = 1
    end
    arrayout = arrayin - mean(arrayin,dim,'omitmissing') ./ std(arrayin,[],dim,'omitmissing');
end