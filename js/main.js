// JavaScript Document
/// <reference path="http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js" />
/// <reference path="numeric-1.2.6.min.js" />



$(document).ready(function () {

	initEvents($("#myCanvas"));

});

// 3次スプライン補間
// 参考：高橋大輔，数値計算，岩波書店，1996.
// 入力 N+1個の点列 points = [[x0,y0],[x1,x2],...,[xn,yn]]
// 出力 各データ点間の区間における多項式の係数を含む2次元配列
// x：x座標値で昇順ソートされた点列のx座標
// y：x座標値で昇順ソートされた点列のy座標
// a：区間[x_i+1 - x_i]における3次式の係数 配列サイズN-1
// b：区間[x_i+1 - x_i]における3次式の係数 配列サイズN-1
// c：区間[x_i+1 - x_i]における3次式の係数 配列サイズN-1
// d：区間[x_i+1 - x_i]における3次式の係数 配列サイズN-1
// ※三次式 y = a*x^3 + b*x^2 + c*x +d
function spline3Interpolation(points) {

	// 入力点数が3未満ではスプライン曲線を定義できないので終了
	if(points.length<3) {
		return null;
	}

	// x座標の値で昇順にソート
	var tmp=numeric.clone(points);
	var sortedPoints=tmp.sort(
		function (a, b) {
			if(a[0]<b[0]) {
				return -1;
			} else {
				return 1;
			}
		}
	);

	// データ点
	var x=new Array(sortedPoints.length);
	var y=new Array(sortedPoints.length);

	// 曲線が定義される区間の数
	var N=points.length-1;

	// 連立一次方程式を構成する行列
	var h=new Array(N);
	var H=numeric.rep([N-1, N-1], 0);
	var u=new Array(N+1);
	var v=new Array(N-1);

	// 各区間の係数 y = a*x^3 + b*x^2 + c*x +d
	var a=new Array(N);
	var b=new Array(N);
	var c=new Array(N);
	var d=new Array(N);

	for(var i=0; i<sortedPoints.length; ++i) {
		x[i]=sortedPoints[i][0];
		y[i]=sortedPoints[i][1];
	}

	for(var i=0; i<N; ++i) {
		h[i] = x[i+1] - x[i];
	}

	H[0][0]=2*(h[0]+h[1]);
	H[0][1]=h[1];
	for(var i=1; i<N-2; ++i) {
		H[i][i-1]=h[i];
		H[i][i]=2*(h[i]+h[i+1]);
		H[i][i+1]=h[i+1];
	}
	H[N-2][N-3]=h[N-2];
	H[N-2][N-2]=2*(h[N-2]+h[N-1]);

	for(var i=1; i<N; ++i) {
		v[i-1]=6*((y[i+1]-y[i])/h[i]-(y[i]-y[i-1])/h[i-1]);
	}

	var uBlock=numeric.solve(H, v);	// from u[1] to u[N-1]
	u[0]=0;
	u[N]=0;
	for(var i=0; i<N-1; ++i) {
		u[i+1] = uBlock[i];
	}

	for(var i=0; i<N; ++i) {
		a[i]=(u[i+1]-u[i])/(6*(x[i+1]-x[i]));
	}

	for(var i=0; i<N; ++i) {
		b[i]=0.5*u[i];
	}

	for(var i=0; i<N; ++i) {
		c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(x[i+1]-x[i])*(2*u[i]+u[i+1])/6;
	}

	for(var i=0; i<N; ++i) {
		d[i]=y[i];
	}

	return [x,y,a,b,c,d];
}

// ラグランジュ補間（と等価な多項式補間）
// 入力 N個の点列 points = [[x0,y0],[x1,x2],...,[xn-1,yn-1]]
// 出力 多項式の係数配列c
//     y = c[0]*x^N-1 + c[1]*x^N-2 + ... + c[N-2]*x + c[N-1]
function lagrangeInterpolation(points) {
	var dim=points.length;
	var A=numeric.rep([dim, dim], 0);
	var b=numeric.rep([dim], 0);
	for(var i=0; i<dim; ++i) {
		for(var j=0; j<dim; ++j) {
			A[i][j]=Math.pow(points[i][0], dim-1-j);
		}
		b[i]=points[i][1];
	}
	var c=numeric.solve(A, b)
	return c;
}


function drawPoints(canvas, points) {
	var context=canvas.get(0).getContext("2d");
	var canvasWidth=canvas.width();
	var canvasHeight=canvas.height();
	context.setTransform(1, 0, 0, 1, 0, 0);
	for(var i=0; i<points.length; ++i) {
		context.beginPath();
		context.arc(points[i][0], points[i][1], 5, 0, 2*Math.PI, true);
		context.fill();
	}
}

function drawSpline(canvas, resultOfSpline) {

	if(resultOfSpline==null) {
		return;
	}

	var px=resultOfSpline[0];
	var py=resultOfSpline[1];
	var a=resultOfSpline[2];
	var b=resultOfSpline[3];
	var c=resultOfSpline[4];
	var d=resultOfSpline[5];

	var context=canvas.get(0).getContext("2d");
	var canvasWidth=canvas.width();
	var canvasHeight=canvas.height();
	var dx=5;
	var numSegments=Math.floor(canvasWidth/dx)+1;
	context.beginPath();
	var x;
	function f(val, id) {
		return a[id]*Math.pow(val-px[id], 3)+b[id]*Math.pow(val-px[id], 2)+c[id]*(val-px[id])+d[id];
	}

	var rangeId=0
	x=0;
	context.moveTo(x, f(x, 0));
	for(var i=1; i<numSegments; ++i) {
		x+=dx;
		if(x>px[rangeId+1]&&rangeId<a.length-1) {
			++rangeId;
		}
		context.lineTo(x, f(x, rangeId));
	}
	context.stroke();
}

function drawPolynomial(canvas, c) {
	var context=canvas.get(0).getContext("2d");
	var canvasWidth=canvas.width();
	var canvasHeight=canvas.height();
	var dx=5;
	var numSegments=Math.floor(canvasWidth/dx)+1;
	context.beginPath();
	var x;
	function f(val) {
		var tmp=0;
		for(var i=0; i<c.length; ++i) {
			tmp+=c[i]*Math.pow(val,c.length-1-i);
		}
		return tmp;
	}
	x=0;
	context.moveTo(x, f(x));
	for(var i=1; i<numSegments; ++i) {
		x+=dx;
		context.lineTo(x, f(x));
	}
	context.stroke();
}

function initEvents(canvas) {

	var canvasWidth=canvas.width();
	var canvasHeight=canvas.height();
	var initPoints = [[200,200], [400,400], [600,200], [800,400]];
	var points=numeric.clone(initPoints);
	var selectPoint=null;

	draw();


	canvas.dblclick(function (event) {
		var canvasOffset=canvas.offset();
		var canvasX=Math.floor(event.pageX-canvasOffset.left);
		var canvasY=Math.floor(event.pageY-canvasOffset.top);
		if(canvasX<0||canvasX>canvasWidth) {
			return;
		}
		if(canvasY<0||canvasY>canvasHeight) {
			return;
		}
		for(var i=0; i<points.length; ++i) {
			if(points[i][0]==canvasX) {
				return;
			}
		}
		points.push([canvasX, canvasY]);
		draw();
	});


	// mouseクリック時のイベントコールバック設定
	canvas.mousedown(function (event) {
		// 左クリック
		if(event.button==0) {
			var canvasOffset=canvas.offset();
			var canvasX=Math.floor(event.pageX-canvasOffset.left);
			var canvasY=Math.floor(event.pageY-canvasOffset.top);
			if(canvasX<0||canvasX>canvasWidth) {
				return;
			}
			if(canvasY<0||canvasY>canvasHeight) {
				return;
			}
			var clickPos=[canvasX, canvasY];
			var dist;
			for(var i=0; i<points.length; ++i) {
				dist=numeric.norm2(numeric.sub(points[i], clickPos));
				if(dist<20) {
					selectPoint=i;
					break;
				}
			}
		} 
	});

	// mouse移動時のイベントコールバック設定
	canvas.mousemove(function (event) {
		var canvasOffset=canvas.offset();
		var canvasX=Math.floor(event.pageX-canvasOffset.left);
		var canvasY=Math.floor(event.pageY-canvasOffset.top);
		if(canvasX<0||canvasX>canvasWidth) {
			return;
		}
		if(canvasY<0||canvasY>canvasHeight) {
			return;
		}
		if(selectPoint!=null) {
			points[selectPoint]=[canvasX, canvasY];
			draw();
		}
	});

	// mouseクリック解除時のイベントコールバック設定
	canvas.mouseup(function (event) {
		selectPoint=null;
		draw();
	});

	// リセットボタン
	$("#reset").click(function () {
		points=numeric.clone(initPoints);
		draw();
	})

	// インプットタグで変化があれば描画更新
	$('input').change(function () {
		draw();
	});

	function draw() {
		var context=canvas.get(0).getContext("2d");
		context.clearRect(0, 0, canvasWidth, canvasHeight);
		drawPoints(canvas, points);
		colors=['red', 'green', 'blue', 'orange', 'cyan'];
		var c;

		context.lineWidth=2;


		if($('#LagrangeCheckBox').is(':checked')) {
			c=lagrangeInterpolation(points);
			context.strokeStyle=colors[0];
			drawPolynomial(canvas, c);
		}

		if($('#SplineCheckBox').is(':checked')) {
			context.strokeStyle=colors[1];
			c=spline3Interpolation(points);
			context.strokeStyle=colors[1];
			drawSpline(canvas, c);
		}
	}


}