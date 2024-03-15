
# this was taken from AutoDock. Use modulo instead??
# Normalizes any radians such that ϑ ∈ [-π, +π]
function normalize_angle(ϑ::Number)
    if ϑ > 3π
        n = (ϑ - π) / (2π)
        ϑ -= 2π*ceil(n)
        normalize_angle(ϑ)
    elseif ϑ < -3π
        n = (-ϑ - π) / (2π)
        ϑ += 2π*ceil(n)
        normalize_angle(ϑ)
    elseif ϑ > π
        ϑ -= 2π
    elseif ϑ < -π
        ϑ += 2π
    end
    ϑ
end


function is_normalized(ϑ::Number)
    ϑ ∈ Interval(π, π) ? true : false
end
