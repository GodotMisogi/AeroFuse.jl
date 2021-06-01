function density(alt)
	if alt <= 36125
		T = 59 - 3.56e-3 * alt
		p = 2116 * ((T + 459.7) / 518.6)^5.256
	elseif 36152 < alt < 82345
		T = -70
		p = 473.1 * exp(1.73 - 4.8e-5 * alt)
	else
		T = -205 + 0.00164 * alt
		p = 51.97 * ((T + 459.7) / 389.98)^(-11.388)
	end
	
	ρ = p / (1718 * (T + 459.7))
end	
