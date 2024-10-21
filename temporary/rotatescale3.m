function [w] = rotatescale3(v0,u)


% %% for 2d
% 
% v0 = [1,2;
%      2,6;
%      4,7;
%      1,9;
%      6,-1];
% 
% v = v0-v0(1,:);
% 
% dx = v(end,1)-v(1,1);
% dy = v(end,2)-v(1,2);
% R = sqrt(dx^2 + dy^2);
% 
% M = [dx, dy;
%      -dy, dx] ./ R^2;
% 
% figure; hold on; xlim([-10 10]); ylim([-10 10]); axis square
% plot(v0(:,1),v0(:,2))
% 
% w = [M*v']'
% 
% plot(w(:,1),w(:,2))

%% for 3d

v = v0-v0(1,:);

[dx,dy,dz] = getRange(v);
R = sqrt(dx^2 + dy^2 + dz^2);

v = v./R;

[dx,dy,~] = getRange(v);

Rxy = sqrt(dx^2 + dy^2);
Mz = [dx, dy, 0;
      -dy, dx, 0;
      0, 0, Rxy] ./ Rxy;
M = Mz;
w = [M*v']';

[dx,~,dz] = getRange(w);
Rxz = sqrt(dx^2 + dz^2);
My = [dx, 0, dz;
      0, Rxz, 0;
      -dz, 0 dx] ./ Rxz;
M = My;
w = [M*w']';

% figure; 
% plot3(v0(:,1),v0(:,2),v0(:,3))
% hold on; xlim([-10 10]); ylim([-10 10]); zlim([-10 10]); axis square
% 
% plot3(w(:,1),w(:,2),w(:,3))s
% 
% xlabel('x')
% ylabel('y')
% zlabel('z')

if exist("u", "var") && ~isempty(u)
    phi = linspace(0,2*pi,1e4)';
    sin_phi = sin(phi);
    cos_phi = cos(phi);
    d = nan(size(phi));
    [dx,dy,dz] = getRange(w);
    R = sqrt(dx^2 + dy^2 + dz^2);
    for i_angle = 1:numel(phi)
        M = [R, 0, 0;
          0, cos_phi(i_angle), -sin_phi(i_angle);
          0, sin_phi(i_angle), cos_phi(i_angle)] ./ R;
        tmp = [M*w']';
        thisd = nanmean(dist(tmp,u));
        d(i_angle) = thisd;
    end
    [~,idx] = nanmin(d); idx = idx(1);
    M = [R, 0, 0;
      0, cos_phi(idx), -sin_phi(idx);
      0, sin_phi(idx), cos_phi(idx)] ./ R;
    w = [M*w']';
end


end

function [d] = dist(v,u)

tmp = v-u;
d = sqrt(sum(tmp.^2));

end

function [dx,dy,dz] = getRange(v)
dx = v(end,1)-v(1,1);
dy = v(end,2)-v(1,2);
dz = v(end,3)-v(1,3);
end












