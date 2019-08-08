function d_thresh=Ransuc_L(temp,tempD)
sampleSize = 3; % number of points to sample per trial
maxDistance = 0.05; % max allowable distance for inliers
points=temp(:,2:3);
fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
evalLineFcn = ...   % distance evaluation function
  @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);

[modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, ...
  sampleSize,maxDistance);
modelInliers = polyfit(points(inlierIdx,1),points(inlierIdx,2),1);
inlierPts = points(inlierIdx,:);
% x = [min(inlierPts(:,1)) max(inlierPts(:,1))];
% y = modelInliers(1)*x + modelInliers(2);
x=tempD(:,2);
y = modelInliers(1)*x + modelInliers(2);
% scatter(x,y);
d_thresh=abs(tempD(:,3)-y);
% plot(abs(points(:,2)-y));
end


% plot(x, y, 'g-')
% legend('Noisy points','Least squares fit','Robust fit');
% hold off