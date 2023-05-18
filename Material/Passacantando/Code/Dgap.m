function v = Dgap(z,alfa,beta)

% D-gap function of a bimatrix game

v = reggap(z,alfa) - reggap(z,beta);

end
