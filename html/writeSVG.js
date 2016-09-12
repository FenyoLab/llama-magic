function writeSVG(path,h,w){
	h = parseInt(h);
	w = parseInt(w);
	var nav = navigator.appName;
	if(nav != "Microsoft Internet Explorer")	{
		document.write('<object data="');
		document.write(path);
		document.write('" height="');
		document.write(h);
		document.write('" width="');
		document.write(w);
		document.write('" type="image/svg+xml" name="SVGEmbed"></object>\n');
	}
	else	{
		document.write('<embed src="');
		document.write(path);
		document.write('" height="');
		document.write(h);
		document.write('" width="');
		document.write(w);
		document.write('" type = "image/svg+xml" PLUGINSPAGE="http://www.adobe.com/svg/viewer/install/" name="SVGEmbed"/>');
	}
}
