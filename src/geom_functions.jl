## calculate the area under the curve
function areaunder(x, y)
    area1 = integrate(abs.(x), abs.(y))

    if area1 > .5
        farea = area1 - .5
    else
        farea = 0
    end

    return farea
end

# calculate the maximum deviation
# note that trajectories have to be normed to end at right side

function maxdev(x, y)

  xchange = x[1] - x[end]
  ychange = y[1] - y[end]
  slope = ychange / xchange
  intercept = y[1] - slope * x[1]

  points = hcat(x, y)
  points_straightline = hcat(x, slope .* x .+ intercept)
  d = sqrt.((points[:,1] .- points_straightline[:,1]).^2 .+ (points[:,2] .- points_straightline[:,2]).^2)

  y_diff = points[:, 2] .- points_straightline[:, 2]
  signs = sign.(y_diff)

  pos_devs = findall(signs .> 0)
  
  if isempty(pos_devs)
    return 0, [0]
  else
    pos_d = d[pos_devs]
    abs_devs = pos_d .* cos.(atan(slope))

    max_dev = maximum(abs_devs)
  end

  return max_dev, pos_devs
end

function interpolate_trajectory(x, y, target_length)
    # Assuming x and y are arrays of the same length
    n = length(x)
    xscale = 1:n
    new_scale = range(1, n, target_length)

    # Create linear interpolations
    x_interp = LinearInterpolation(xscale, x)
    y_interp = LinearInterpolation(xscale, y)

    # Interpolate to new scale
    new_x = [x_interp(xi) for xi in new_scale]
    new_y = [y_interp(yi) for yi in new_scale]

    return new_x, new_y
end