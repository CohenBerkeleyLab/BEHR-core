function [ x, y ] = uncross_pix_corners( x, y )
%[X, Y] = UNCROSS_PIX_CORNERS(X, Y)
%   Given vectors of four pixel corners, this function will ensure
%   that the corners do not "cross," that is, points 2 and 3 are not
%   catty-corner from one another, nor are points 1 and 4. It will
%   return uncrossed corners.
%
%   Josh Laughner <joshlaugh5@gmail.com> Mar 2016

m1 = (y(3) - y(2))/(x(3) - x(2));
b1 = y(2) - m1*x(2);
m2 = (y(4) - y(1))/(x(4) - x(1));
b2 = y(1) - m2*x(1);

flip_bool = false;
if ~isinf(m1) && ~isinf(m2)
    % As long as neither slope is infinite, solve for the x-coordinate
    % of the intercept and see if it falls inside the polygon - if so,
    % the corners need flipped.
    inpt = (b2-b1)/(m1-m2);
    if inpt > min(x(2:3)) && inpt < max(x(2:3))
        flip_bool = true;
    end
elseif isinf(m1) && ~isinf(m2)
    % If one is infinite, fine the y-coord where the other one is at
    % it's x-coordinate and do the same test
    inpt = m2*x(2)+b2;
    if inpt >= min(y([1,4])) && inpt <= max(y([1,4]))
        flip_bool = true;
    end
elseif isinf(m2) && ~isinf(m1)
    inpt = m1*x(1) + b1;
    if inpt >= min(y(2:3)) && inpt <= max(y(2:3))
        flip_bool = true;
    end
    % If both are infinite, they are parallel and the corners do not
    % need flipped.
end

if flip_bool
    tmp = x(4);
    x(4) = x(3);
    x(3) = tmp;
    tmp = y(4);
    y(4) = y(3);
    y(3) = tmp;
end

m1 = (y(2) - y(1))/(x(2) - x(1));
b1 = (y(1) - m1*x(1));
m2 = (y(4) - y(3))/(x(4) - x(3));
b2 = (y(3) - m2*x(3));

flip_bool = false;
if ~isinf(m1) && ~isinf(m2)
    % As long as neither slope is infinite, solve for the x-coordinate
    % of the intercept and see if it falls inside the polygon - if so,
    % the corners need flipped.
    inpt = (b2-b1)/(m1-m2);
    if inpt > min(x(1:2)) && inpt < max(x(1:2))
        flip_bool = true;
    end
elseif isinf(m1) && ~isinf(m2)
    % If one is infinite, fine the y-coord where the other one is at
    % it's x-coordinate and do the same test
    inpt = m2*x(2)+b2;
    if inpt >= min(y(3:4)) && inpt <= max(y(3:4))
        flip_bool = true;
    end
elseif isinf(m2) && ~isinf(m1)
    inpt = m1*x(1) + b1;
    if inpt >= min(y(1:2)) && inpt <= max(y(1:2))
        flip_bool = true;
    end
    % If both are infinite, they are parallel and the corners do not
    % need flipped.
end

if flip_bool
    tmp = x(2);
    x(2) = x(3);
    x(3) = tmp;
    tmp = y(2);
    y(2) = y(3);
    y(3) = tmp;
end

end


