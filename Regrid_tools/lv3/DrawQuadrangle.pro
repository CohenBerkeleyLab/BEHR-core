function calcline, y,x1,y1,x2,y2
	if y2 eq y1 then return, x1 else return, x1+round(DOUBLE(y-y1)*DOUBLE(x2-x1)/DOUBLE(y2-y1))
end

pro clip, value, minvalue, maxvalue
	if value lt minvalue then value=minvalue
	if value ge maxvalue then value=maxvalue-1

end

pro Quadrangle, array, coord, value
;written by Mark Wenig
;fills all cells of 'array' between 'coord' with 'value'
	x1=round(coord(0,0))
	y1=round(coord(1,0))
	x2=round(coord(0,1))
	y2=round(coord(1,1))
	x3=round(coord(0,2))
	y3=round(coord(1,2))
	x4=round(coord(0,3))
	y4=round(coord(1,3))
	if y2 lt y1 then ExchangeCoord, x1,y1,x2,y2
	if y3 lt y1 then ExchangeCoord, x1,y1,x3,y3
	if y4 lt y1 then ExchangeCoord, x1,y1,x4,y4
	if y2 gt y3 then ExchangeCoord, x2,y2,x3,y3
	if y4 gt y3 then ExchangeCoord, x4,y4,x3,y3
	if x4 gt x2 then ExchangeCoord, x4,y4,x2,y2
	bottom=y1
	top=y3-1
	dim = SIZE( array , /DIMENSIONS)
	if (bottom lt dim(1)) and (top ge 0) then begin
		clip, bottom, 0, dim(1)
		clip, top, 0, dim(1)

		for y = bottom, top,1 do begin
			if (x2 ge calcline(y2,x1,y1,x3,y3)) and (x4 le calcline(y4,x1,y1,x3,y3)) then begin
				if y lt y4 then left=calcline(y,x1,y1,x4,y4)$
				else left=calcline(y,x4,y4,x3,y3)
				if y lt y2 then right=calcline(y,x1,y1,x2,y2)$
				else right=calcline(y,x2,y2,x3,y3)
			end else begin
				left=calcline(y,x1,y1,x3,y3)
				if y2 gt y4 then ExchangeCoord, x4,y4,x2,y2
				if y lt y2 then right=calcline(y,x1,y1,x2,y2)$
				else if y lt y4 then right=calcline(y,x2,y2,x4,y4)$
				else right=calcline(y,x4,y4,x3,y3)
				if (x2 le calcline(y2,x1,y1,x3,y3)) and (x4 le calcline(y4,x1,y1,x3,y3)) then begin
					dummy=left
					left=right
					right=dummy
				end
			end
			right-=1
			if (right ge 0) and (left lt dim(0)) then begin
				clip, left, 0, dim(0)
				clip, right, 0, dim(0)
				for x = left, right do array(x,y)+=value
			end
		end
	end
end
