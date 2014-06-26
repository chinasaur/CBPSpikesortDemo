function bool = haveipt()
v = ver();
bool = any(strcmp('Image Processing Toolbox', {v.Name}));