function warning_assert(expr, msg)
	if ~expr
		warning(msg);
	end
end