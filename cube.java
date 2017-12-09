import java.awt.Canvas;
import javax.swing.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Insets;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.awt.image.BufferedImage;

public class cube extends JPanel {

    int maxX, maxY, centerX, centerY, edgelength, imagesize;
    float rwidth = 100, rheight = 100, pixelsizex, pixelsizey, pixelsize;
    float rho, theta, phi;
    float sint, sinp, cost, cosp;
    float m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44;
    float ex, ey, ez, objectsize = 12,objectsize1=12;
    float sx, sy;
    float d1,d2;
    int im = 1;
    float a0, a1, a2, a3, a4, a5, a6, a7;
    Point2D[] sc,sc1;
    float xe = 4, ye = 10, ze = 2;
    float xk = 4, yk = 10, zk = 2;
    Image image,image1;
    Graphics g2,g3;
    int w, h;
    Point3D mvc = new Point3D(d1, d1, d1);
    Point3D[] wc, e,e1;
    float p1x, p1y, p1z, p2x, p2y, p2z;
    double iscrx, iscry, jscrx, jscry;
    int one, two, three;
    Point2D pscr;
    Point2D qscr;
    int istart;
    int k;
    Triangle[] tr = new Triangle[12];
    Triangle[] all = new Triangle[12];
    int[] colorcode = new int[12];
    Triangle[] tri = new Triangle[12];
    float za1, zb1, zc1, zaverage;
    float[] ztr = new float[12];
    float p1x1, p1y1, p1z1, p2x1, p2y1, p2z1;
    double sunx = -1 / Math.sqrt(3), suny = 1 / Math.sqrt(3), sunz = 1 / Math.sqrt(3);
    double inprod, inprodMin = 1e30, inprodMax = -1e30, inprodRange;
    Point2D[] pointsa, pointsb, pointsc, vscr;
    float cax, cay, cbx, cby, ccx, ccy;    
    int x1, y1;
    int cx, cy;
    float xx, yy, zz;
    int ki = 1;
    int[] color = new int[12];
    private BufferedImage bf;
    
    void Dime() {
        Dimension d = getSize();
        maxX = d.width - 1;
        maxY = d.height - 1;
        if (im == 1) {
            centerX = (int) (maxX / (3.5));
            centerY = (int) (maxY / 1.64);   //1.47,1.60
            im++;
            cx = centerX;
            cy = centerY;
        }
        imagesize = Math.min(maxX / 4, maxY / 4);

    }

    public cube() {
        

        addMouseListener(new MouseAdapter() {

            public void mouseClicked(MouseEvent e) {
                Dimension d = getSize();
                maxX = d.width - 1;
                maxY = d.height - 1;

                x1 = centerX;
                y1 = centerY;
                animate();
            }
        });
    }
  
    void animate() {
        Thread thread = new Thread() {
            @Override
            public void run() {
                try {
                    while (x1 <= (int) (maxX / 1.47)) {
                        x1 = (int) (x1 + 1);
                        float x11 = fx(x1);
                        float y11 = fy(y1);
                        System.out.println("maxX" + maxX);
                        System.out.println("maxY" + maxY);
                        tranx(x11, y11);
                        centerX = x1;
                        centerY = y1;
                        repaint();
                        Thread.sleep(10);
                    }

                } catch (InterruptedException ex) {
                }
            }
        };
        thread.start();
    }

        
    int ix(float x) { return Math.round(centerX + x);}

    int iy(float y) {return Math.round(centerY - y);}

    float fx(int X) { return X - centerX;}

    float fy(int Y) {return centerY - Y;}
    void tranx(float X, float Y) {

        float xd = X;
        float yd = Y;
        float z3d = -rho;
        float x3d = -((xd * z3d) / (d1));
        float y3d = -((yd * z3d) / (d1));

        float b11 = (float) -Math.sin(theta), b12 = (float) Math.cos(theta), b13 = 0, b14 = 0,
                //b21 = (float) (-Math.cos(theta)*Math.sin(phi)),
                b21 = (float) (-Math.cos(theta) * Math.cos(phi)),
                b22 = (float) (-Math.sin(theta) * Math.cos(phi)),
                b23 = (float) (Math.sin(phi)), b24 = 0,
                //b31 = (float) (Math.cos(phi)*Math.sin(phi)),
                b31 = (float) (Math.cos(theta) * Math.sin(phi)),
                b32 = (float) (Math.sin(theta) * Math.sin(phi)),
                b33 = (float) (Math.cos(phi)),
                b34 = 0,
                //b41 = (float) (Math.cos(theta)*Math.sin(theta)*rho),
                b41 = (float) (Math.cos(theta) * Math.sin(phi) * rho),
                b42 = (float) (Math.sin(theta) * Math.sin(phi) * rho),
                b43 = (float) (Math.cos(phi) * rho),
                b44 = 1;
        float x2 = ((x3d * b11) + (y3d * b21) + (z3d * b31) + (b41));
        float y2 = ((x3d * b12) + (y3d * b22) + (z3d * b32) + (b42));
        float z2 = ((x3d * b13) + (y3d * b23) + (z3d * b33) + (b43));
        float ww = 1;
   
        xx = x2;
        yy = y2;
        zz = z2;
        xe = xe - xx;
        ye = ye - yy;
        ze = ze - zz;
        rho = (float) Math.sqrt(((xe) * (xe)) + ((ye) * (ye)) + ((ze) * (ze)));

        theta = (float) Math.atan2(ye, xe);
        phi = (float) Math.acos(ze / rho);
        System.out.println("xx: " + xx);
        System.out.println("yy: " + yy);
        System.out.println("zz: " + zz);
        System.out.println("xe: " + xe);
        System.out.println("ye: " + ye);
        System.out.println("ze: " + ze);

    }


  public void update(Graphics g)
	{
		paint(g);
	}

   
    public void paint(Graphics g) {

        e = new Point3D[8];
        Dime();
        Dimension d = getSize();

        if (w != d.width || h != d.height) {
            w = d.width;
            h = d.height;
            image = createImage(w, h);
            image1=createImage(w,h);
            g2 = image.getGraphics();
            g3=image1.getGraphics();
            
        }
        g2.clearRect(0, 0, w, h);
        g3.clearRect(0,0,w,h);
        dcube();
        Point2D[] f1, f2, f3, f4, f5, f6;
        f1 = new Point2D[4];
        f1[0] = sc[0];
        f1[1] = sc[1];
        f1[2] = sc[5];
        f1[3] = sc[4];
        triangulate(f1, tr);
        all[0] = tr[0];
        all[1] = tr[1];
        f2 = new Point2D[4];
        f3 = new Point2D[4];
        f4 = new Point2D[4];
        f5 = new Point2D[4];
        f6 = new Point2D[4];
        f2[0] = sc[1];
        f2[1] = sc[2];
        f2[2] = sc[6];
        f2[3] = sc[5];
        triangulate(f2, tr);
        all[2] = tr[0];
        all[3] = tr[1];
        f3[0] = sc[0];
        f3[1] = sc[3];
        f3[2] = sc[7];
        f3[3] = sc[4];
        triangulate(f3, tr);
        all[4] = tr[0];
        all[5] = tr[1];
        f5[0] = sc[0];
        f5[1] = sc[1];
        f5[2] = sc[2];
        f5[3] = sc[3];
        triangulate(f5, tr);
        all[6] = tr[0];
        all[7] = tr[1];
        f6[0] = sc[3];
        f6[1] = sc[2];
        f6[2] = sc[6];
        f6[3] = sc[7];
        triangulate(f6, tr);
        all[8] = tr[0];
        all[9] = tr[1];

        for (int i = 0; i < 8; i++) {

            line(g2, 2, 3,sc);
            line(g2, 0, 1,sc);
            line(g2, 1, 5,sc);
            line(g2, 5, 4,sc);
            line(g2, 4, 0,sc);
            line(g2, 7, 3,sc);
            line(g2, 7, 6,sc);
            line(g2, 6, 5,sc);
            line(g2, 0, 3,sc);
            line(g2, 1, 2,sc);
            line(g2, 4, 7,sc);
            line(g2, 2, 6,sc);
            
            
            line(g2, 0,1,sc1);
            line(g2, 0,2,sc1);
            line(g2, 2,3,sc1);
            line(g2, 3,1,sc1);
            line(g2, 2,4,sc1);
            line(g2, 3,5,sc1);
            line(g2,1,6,sc1);
            line(g2,4,5,sc1);
            line(g2,5,6,sc1);
            g.drawImage(image, 0, 0, null);
        }
            

        tri[0] = new Triangle(e[0], e[1], e[5]);
        tri[1] = new Triangle(e[5], e[4], e[0]);
        tri[2] = new Triangle(e[1], e[2], e[6]);
        tri[3] = new Triangle(e[6], e[5], e[1]);
        tri[4] = new Triangle(e[0], e[3], e[7]);
        tri[5] = new Triangle(e[7], e[4], e[0]);
        tri[6] = new Triangle(e[4], e[5], e[6]);
        tri[7] = new Triangle(e[6], e[7], e[4]);
        tri[8] = new Triangle(e[0], e[1], e[2]);
        tri[9] = new Triangle(e[2], e[3], e[0]);
        tri[10] = new Triangle(e[3], e[2], e[6]);
        tri[11] = new Triangle(e[6], e[7], e[3]);

        
        for (int t = 0; t < 12; t++) {
            za1 = tri[t].a3.z;
            zb1 = tri[t].b3.z;
            zc1 = tri[t].c3.z;
            zaverage = za1 + zb1 + zc1;
            ztr[t] = zaverage;
            p1x1 = (tri[t].b3.x) - (tri[t].a3.x);
            p1y1 = (tri[t].b3.y) - (tri[t].a3.y);
            p1z1 = (tri[t].b3.z) - (tri[t].a3.z);
            p2x1 = (tri[t].c3.x) - (tri[t].b3.x);
            p2y1 = (tri[t].c3.y) - (tri[t].b3.y);
            p2z1 = (tri[t].c3.z) - (tri[t].b3.z);
            
            double a = ((p1y1 * p2z1) - (p1z1 * p2y1)), b = ((p1z1 * p2x1) - (p1x1 * p2z1)), c = ((p1x1 * p2y1) - (p1y1 * p2x1)),
                    h = (a * tri[t].c3.x) + (b * tri[t].c3.y) + (c * tri[t].c3.z);
            
            int ccode = colorcode(a, b, c);
		        int de;
		        if((t>5 && t<8))
		        	de = 2;
		        else
		        	de = 1;
		        colorcode[t] = ccode;
		        System.out.println("colorcode: "+colorcode[t]+ ", triangle: " +t);
		        color[t] = de;
        }
         sort(tri, colorcode, ztr, 0, 11, color);
			 for (int iTria=0; iTria<12; iTria++)
		      {  
				 for(int yy=0; yy<8; yy++)
				 {
					 if(((-d1*(tri[iTria].a3.x/tri[iTria].a3.z)) == sc[yy].x)&&
							 ((-d1*(tri[iTria].a3.y/tri[iTria].a3.z)) == sc[yy].y))
					 {
						
						cax = sc[yy].x;
						cay = sc[yy].y;
					 }
				 }
				 for(int yy=0; yy<8; yy++)
				 {
					 if(((-d1*(tri[iTria].b3.x/tri[iTria].b3.z)) == sc[yy].x)&&
							 ((-d1*(tri[iTria].b3.y/tri[iTria].b3.z)) == sc[yy].y))
					 {
						 
						 cbx = sc[yy].x;
						 cby = sc[yy].y;
					 }
				 }
				for(int yy=0; yy<8; yy++)
				{
					 if(((-d1*(tri[iTria].c3.x/tri[iTria].c3.z)) == sc[yy].x)&&
							 ((-d1*(tri[iTria].c3.y/tri[iTria].c3.z)) == sc[yy].y))
					 {
						 
						 ccx = sc[yy].x;
						 ccy = sc[yy].y;
					 }
				}
				
		         int ccode = colorcode[iTria];
		         //System.out.println("ccode: " +ccode);
		         int de = color[iTria];
		         if(de == 1)
		         {
		        	 //System.out.println("blue");
		        	 //g.setColor(Color.BLUE);
		        	 g.setColor(new Color(0,0,ccode));
		         }
		        if(de == 2)
		        {
		        	 //System.out.println("red");
		        	 //g.setColor(Color.RED);
		        	g.setColor(new Color(ccode,0,0));
                        }
		         int[] x = {ix(cax), ix(cbx), ix(ccx)};
		         int[] y = {iy(cay), iy(cby), iy(ccy)};
		         g.fillPolygon(x, y, 3);
		         
		      }

    }
    
    public void line(Graphics g2, int k1, int k2,Point2D[] sc3) {

        Point3D p;
        p = e[k1];
        Point3D q;
        q = e[k2];

        float px, py, qx, qy;
        px = sc3[k1].x;
        py = sc3[k1].y;
        qx = sc3[k2].x;
        qy = sc3[k2].y;
        Point2D pscr = new Point2D(px, py);
        Point2D qscr = new Point2D(qx, qy);

        linesegment(g2, p, q, pscr, qscr, k1, k2, istart,sc3);
    }
    
    public void linesegment(Graphics g2, Point3D p, Point3D q, Point2D pscr, Point2D qscr, int iP, int iQ, int istart,Point2D[] sc) {
        // TODO Auto-generated method stub
        //System.out.println("pscr, sc[k1.x]" +pscr.x);

        double u1 = qscr.x - pscr.x, u2 = qscr.y - pscr.y;
        double minx = Math.min(pscr.x, qscr.x);
        double maxx = Math.max(pscr.x, qscr.x);
        double miny = Math.min(pscr.y, qscr.y);
        double maxy = Math.min(pscr.y, qscr.y);
        double zP = p.z, zQ = q.z;
        double minz = Math.min(zP, zQ);
        for (int t = istart; t < 10; t++) {
            float a1, a2, b1, b2, c1, c2;
            a1 = all[t].a.x;
            a2 = all[t].a.y;
            b1 = all[t].b.x;
            b2 = all[t].b.y;
            c1 = all[t].c.x;
            c2 = all[t].c.y;
            Point2D ascr = new Point2D(a1, a2);
            Point2D bscr = new Point2D(b1, b2);
            Point2D cscr = new Point2D(c1, c2);

            for (int o = 0; o < 8; o++) {
                if (all[t].a.x == sc[o].x && all[t].a.y == sc[o].y) {
                    one = o;
                }
            }
            for (int o = 0; o < 8; o++) {
                if (all[t].b.x == sc[o].x && all[t].b.y == sc[o].y) {
                    two = o;
                }

            }
            for (int o = 0; o < 8; o++) {
                if (all[t].c.x == sc[o].x && all[t].c.y == sc[o].y) {
                    three = o;
                }
            }

            double eps = 0.1; // Relative to numbers of pixels
            if (area2(ascr, bscr, pscr) < eps
                    && area2(ascr, bscr, qscr) < eps
                    || area2(bscr, cscr, pscr) < eps
                    && area2(bscr, cscr, qscr) < eps
                    || area2(cscr, ascr, pscr) < eps
                    && area2(cscr, ascr, qscr) < eps) {
                //System.out.println("test 4");
                continue;
            }

            //test 5
            double pqa = area2(pscr, qscr, ascr);
            double pqb = area2(pscr, qscr, bscr);
            double pqc = area2(pscr, qscr, cscr);

            if (pqa < +eps && pqb < +eps && pqc < +eps
                    || pqa > -eps && pqb > -eps && pqc > -eps) {
                //System.out.println("test 5");
                continue;
            }

            //test 6
            p1x = (e[two].x) - (e[one].x);
            p1y = (e[two].y) - (e[one].y);
            p1z = (e[two].z) - (e[one].z);
            p2x = (e[three].x) - (e[two].x);
            p2y = (e[three].y) - (e[two].y);
            p2z = (e[three].z) - (e[two].z);
            //System.out.println("p1x, p1y, p1z: " +p1x+ ", " +p1y+ " , " +p1z);
            //System.out.println("p2x, p2y, p2z: " +p2x+ ", " +p2y+ " , " +p2z);
            double a = ((p1y * p2z) - (p1z * p2y)), b = ((p1z * p2x) - (p1x * p2z)), c = ((p1x * p2y) - (p1y * p2x)),
                    h = (a * e[three].x) + (b * e[three].y) + (c * e[three].z), eps1 = 1e-5 * Math.abs(h),
                    hP = a * p.x + b * p.y + c * p.z,
                    hQ = a * q.x + b * q.y + c * q.z;
            if (hP > h - eps1 && hQ > h - eps1) {
                //System.out.println("test 6");
                continue;
            }

            //test 7
            boolean pInside = insideTriangle(ascr, bscr, cscr, pscr);
            boolean qInside = insideTriangle(ascr, bscr, cscr, qscr);
            if (pInside && qInside) {
                //System.out.println("test 7");
                return;
            }

            //test 8
            double h1 = h + eps1;
            if (hP > h1 && pInside || hQ > h1 && qInside) {
                //System.out.println("test 8");
                continue;
            }

            // test 9
            double lambdaMin = 1.0, lambdaMax = 0.0;
            for (int ii = 0; ii < 3; ii++) {
                double v1 = bscr.x - ascr.x, v2 = bscr.y - ascr.y,
                        w1 = ascr.x - pscr.x, w2 = ascr.y - pscr.y,
                        denom = u2 * v1 - u1 * v2;
                if (denom != 0) {
                    double mu = (u1 * w2 - u2 * w1) / denom;

                    if (mu > -0.0001 && mu < 1.0001) {
                        double lambda = (v1 * w2 - v2 * w1) / denom;

                        if (lambda > -0.0001 && lambda < 1.0001) {
                            if (pInside != qInside
                                    && lambda > 0.0001 && lambda < 0.9999) {
                                lambdaMin = lambdaMax = lambda;
                                break;

                            }
                            if (lambda < lambdaMin) {
                                lambdaMin = lambda;
                            }
                            if (lambda > lambdaMax) {
                                lambdaMax = lambda;
                            }
                        }
                    }
                }
                Point2D temp = ascr;
                ascr = bscr;
                bscr = cscr;
                cscr = temp;

            }

            if (!pInside && lambdaMin > 0.001) {
                double iScrx = pscr.x + lambdaMin * u1,
                        iScry = pscr.y + lambdaMin * u2;
                double zI = 1 / (lambdaMin / zQ + (1 - lambdaMin) / zP),
                        xI = -zI * iScrx / d1, yI = -zI * iScry / d1;

                if (a * xI + b * yI + c * zI > h1) {
                    continue; // This triangle does not obscure PQ.
                }
                Point2D iScr = new Point2D((float) iScrx, (float) iScry);

                if (distance2(iScr, pscr) >= 1.0) {
                    linesegment(g2, p, new Point3D((float) xI, (float) yI, (float) zI), pscr, iScr, iP, -1, t + 1,sc);
                }

            }
            if (!qInside && lambdaMax < 0.999) {
                double jscrx = pscr.x + lambdaMax * u1,
                        jscry = pscr.y + lambdaMax * u2;
                double zJ
                        = 1 / (lambdaMax / zQ + (1 - lambdaMax) / zP),
                        xJ = -zJ * jscrx / d1, yJ = -zJ * jscry / d1;
                if (a * xJ + b * yJ + c * zJ > h1) {
                    continue;
                }
                Point2D jScr = new Point2D((float) jscrx, (float) jscry);
                if (distance2(jScr, qscr) >= 1.0) {
                    linesegment(g2, q, new Point3D((float) xJ, (float) yJ, (float) zJ), qscr, jScr, iQ, -1, t + 1,sc);
                }
            }
            return;

        }
        //System.out.println("helloo");
        drawLine(g2, pscr.x, pscr.y, qscr.x, qscr.y);

    }

    private void drawLine(Graphics g2, float x1, float y1, float x2, float y2) {
        // TODO Auto-generated method stub
        if (x1 != x2 || y1 != y2) {
            //g2.setColor(Color.GREEN);
            g2.drawLine(ix(x1), iy(y1), ix(x2), iy(y2));
           
        }

    }

    private void sort(Triangle[] tri, int[] colorcode, float[] ztr, int l, int r, int[] color) {
		// TODO Auto-generated method stub
		int i = l, j = r,  wInt, wd;
	      float x = ztr[(i + j)/2], w;
	      Triangle wTria;
	      do
	      {  while (ztr[i] < x) i++;
	         while (ztr[j] > x) j--;
	         if (i < j)
	         {  w = ztr[i]; ztr[i] = ztr[j]; ztr[j] = w;
	            wTria = tri[i]; tri[i] = tri[j]; tri[j] = wTria;
	            wInt = colorcode[i]; colorcode[i] = colorcode[j];
	            colorcode[j] = wInt;
	            wd = color[i]; color[i] = color[j];
	            color[j] = wd;
	            i++; j--;
	         } else
	      if (i == j) {i++; j--;}
	   } while (i <= j);
	   if (l < j) sort(tri, colorcode, ztr, l, j, color);
	   if (i < r) sort(tri, colorcode, ztr, i, r, color);
	   for(int kk=0; kk<12; kk++ )
	   {
		   //System.out.println("ztr: " +ztr[kk]+ ", kk: " +kk);
	   }
	}

	private int colorcode(double a, double b, double c) {
		// TODO Auto-generated method stub
		inprod = a * sunx + b * suny + c * sunz;
		System.out.println("inprod: " +inprod);
		System.out.println("inprodMin: " +inprodMin);
		System.out.println("inprodMax: " +inprodMax);
		if (inprod < inprodMin) inprodMin = inprod; 
        if (inprod > inprodMax) inprodMax = inprod;
		inprodRange = inprodMin - inprodMax;
		int bb = (int)Math.round(((inprod - inprodMin)/inprodRange)*255);
		System.out.println("bb: "+bb);
		return Math.abs((int)Math.round(((inprod - inprodMin)/inprodRange)*255));
	}
    
    
    public double distance2(Point2D i11, Point2D j11) {
        // TODO Auto-generated method stub
        {
            float dx = i11.x - j11.x,
                    dy = i11.y - j11.y;
            return dx * dx + dy * dy;
        }
    }

    public void triangulate(Point2D[] f1, Triangle[] tr) {
        // TODO Auto-generated method stub
        int n = 4, j = n - 1, i1 = 0, i2, i3;
        int[] ver = new int[n];

        for (int i = 0; i < n; i++) {
            ver[j] = i;
            j = i;
        }
        for (k = 0; k < n - 2; k++) {  // Find a suitable triangle, consisting of two edges
            // and an internal diagonal:
            Point2D a, b, c;
            boolean triaFound = false;
            int count = 0;
            while (!triaFound && ++count < n) {
                i2 = ver[i1];
                i3 = ver[i2];
                a = f1[i1];
                b = f1[i2];
                c = f1[i3];
                if (area2(a, b, c) >= 0) {  // Edges AB and BC; diagonal AC.
                    // Test to see if no other polygon vertex
                    // lies within triangle ABC:
                    j = ver[i3];
                    while (j != i1 && !insideTriangle(a, b, c, f1[j])) {
                        j = ver[j];
                    }
                    if (j == i1) {  // Triangle ABC contains no other vertex:
                        tr[k] = new Triangle(a, b, c);
                        ver[i1] = i3;
                        triaFound = true;
                    }
                }
                i1 = ver[i1];
            }
            //if (count == n)
            //{  System.out.println("Not a simple polygon" +
            //     " or vertex sequence not counter-clockwise.");
            //   System.exit(1);
            //}
        }
    }

    public boolean insideTriangle(Point2D a, Point2D b, Point2D c, Point2D p) {
        // TODO Auto-generated method stub
        return area2(a, b, p) >= 0
                && area2(b, c, p) >= 0
                && area2(c, a, p) >= 0;

    }

    public float area2(Point2D a, Point2D b, Point2D c) {
        // TODO Auto-generated method stub
        return ((a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x));
    }

    public void dcube() {
        // TODO Auto-generated method stub
        //Point3D[] wc;
        wc = new Point3D[8];
        sc = new Point2D[8];
        sc1 = new Point2D[8];

        wc[0] = new Point3D(mvc.x + 1, mvc.y - 1, mvc.z - 1);
        wc[1] = new Point3D(mvc.x + 1, mvc.y + 1, mvc.z - 1);
        wc[2] = new Point3D(mvc.x - 1, mvc.y + 1, mvc.z - 1);
        wc[3] = new Point3D(mvc.x - 1, mvc.y - 1, mvc.z - 1);
        wc[4] = new Point3D(mvc.x + 1, mvc.y - 1, mvc.z + 1);
        wc[5] = new Point3D(mvc.x + 1, mvc.y + 1, mvc.z + 1);
        wc[6] = new Point3D(mvc.x - 1, mvc.y + 1, mvc.z + 1);
        wc[7] = new Point3D(mvc.x - 1, mvc.y - 1, mvc.z + 1);
        if (ki == 1) {
            rho = (float) Math.sqrt(((xe * xe)) + (ye * ye) + (ze * ze));
            //rho = 15;
            phi = (float) Math.acos(ze / rho);
            theta = (float) Math.atan2(ye, xe);
            ki++;
        }
        objectsize = (float) Math.sqrt(25F);
        //objectsize1 = (float) Math.sqrt(125F);
        //rho = 5*objectsize;
        d1 = rho * (imagesize / objectsize);
        //d2=rho * (imagesize / objectsize1);
        //ze = -d1*(fx(x1)/x1);
        eye(wc, d1);
        //eye(wc,d2);
    }



    class Point2D {

        float x, y;

        Point2D(float x, float y) {
            this.x = x;
            this.y = y;
        }
    }

    public void eye(Point3D[] wc, float d1) {
        // TODO Auto-generated method stub
         Dimension d = getSize();
        maxX = d.width - 1;
        maxY = d.height - 1;
        
        sint = (float) Math.sin(theta);
        cost = (float) Math.cos(theta);
        sinp = (float) Math.sin(phi);
        cosp = (float) Math.cos(phi);
        m11 = -sint;
        m12 = ((-cosp) * (cost));
        m13 = ((sinp) * (cost));
        m14 = 0;
        m21 = cost;
        m22 = ((-cosp) * (sinp));
        m23 = ((sinp) * (sint));
        m24 = 0;
        m31 = 0;
        m32 = sinp;
        m33 = cosp;
        m34 = 0;
        m41 = 0;
        m42 = 0;
        m43 = -rho;
        m44 = 1;
        for (int i = 0; i < 8; i++) {
            Point3D w = wc[i];
            Point3D w1 = wc[i];
            ex = ((w.x * m11) + (w.y * m21) + (w.z * m31) + (m41));
            ey = ((w.x * m12) + (w.y * m22) + (w.z * m32) + (m42));
            ez = ((w.x * m13) + (w.y * m23) + (w.z * m33) + (m43));
            sx = (-d1 * (ex / ez));
            sy = (-d1 * (ey / ez));
            //System.out.println("d1 here: "+d1);
            e[i] = new Point3D(ex, ey, ez);
            //e[i+10]=new Point3D(ex, ey, ez);
            //System.out.println("e[i]: " +e[i].x+ ", " +e[i].y+ ", " +i);
            //Point2D sc = new Point2D();
            //sc[i] = new Point2D(sx, sy);
            sc[i] = new Point2D(sx, sy);

        }
        sc1[0]=new Point2D(fx(maxX/4),fy(maxY/2));
        sc1[1]=new Point2D(fx(3*maxX/4),fy(maxY/2));
        sc1[2]=new Point2D(fx(maxX/6),fy(3*maxY/4));
        sc1[3]=new Point2D(fx(5* maxX / 7),fy (3 * maxY / 4));
        sc1[4]=new Point2D(fx(maxX / 6), fy(8 * maxY / 9));
        sc1[5]=new Point2D(fx(5 * maxX / 7), fy(8 * maxY / 9));
        sc1[6]=new Point2D(fx(3 * maxX / 4), fy(2 * maxY / 3));
        sc1[7]=new Point2D(fx(maxX / 6), fy(8 * maxY / 9));
    }
    
    
    

}
