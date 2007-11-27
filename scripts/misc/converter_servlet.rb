#!/usr/bin/env ruby

#dodgy servlet to interpolate between 2 points

require 'webrick'
require 'cgi'

include WEBrick

def deg(x1,x2,y1,y2,z)
	return y1+(y2-y1)*(z-x1).to_f/(x2-x1).to_f
end

class MapServlet < HTTPServlet::AbstractServlet

    def do_GET(req, resp)

    resp['Content-type']='text/html'
    if req.query['x1']
	    x1=req.query['x1'].to_f
	    x2=req.query['x2'].to_f
	    y1=req.query['y1'].to_f
	    y2=req.query['y2'].to_f
	    z=req.query['z'].to_f

	    resp.body=<<TEMP;
<html>
<form action="physmap" method="get">
<p>Gene A: physical position: <input name="x1" type="text" value="#{x1}"></input>bp genetic position <input name="y1" value="#{y1}"></input>cM </p>
<p>Gene B: physical position: <input name="z" type="text" value="#{z}"></input>bp genetic position #{deg(x1,x2,y1,y2,z)} cM </p>
<p>Gene C: physical position: <input name="x2" type="text" value="#{x2}"></input>bp genetic position <input name="y2" value="#{y2}"></input>cM </p>
<input type="submit">
<input type="button" name="reset" value="reset" onclick="window.location = 'physmap' ">
</form>
</html>
TEMP
      raise HTTPStatus::OK
    else
resp.body=<<TEMP;
<html>
<form action="physmap" method="get">
<p>Gene A: physical position: <input name="x1" type="text"></input>bp genetic position <input name="y1"></input>cM </p>
<p>Gene B: physical position: <input name="z" type="text"></input>bp</p>
<p>Gene C: physical position: <input name="x2" type="text"></input>bp genetic position <input name="y2"></input>cM </p>
<input type="submit">
</form>
</html>
TEMP

      raise HTTPStatus::OK
    end
  end

end

def start_webrick(config={})
	config.update(:Port => 9988)
	server = HTTPServer.new(config)
	yield server if block_given?
	['INT','TERM'].each {|signal|
		trap(signal){server.shutdown}
	}
	server.start
end


start_webrick {|s|
	s.mount('/physmap', MapServlet) 
}

