function xNew = stateEvolution(x)
theta = 30;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

xNew = R*x;

end