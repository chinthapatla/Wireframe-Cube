import java.awt.Frame;
import java.awt.PopupMenu;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;


public class wireframe extends Frame {
	public static void main(String args[])
	{
		new wireframe();
	}
	
	public wireframe() {
		// TODO Auto-generated method stub
		Frame f = new Frame("Final Lab");
        f.setSize(1000, 1000);
        f.addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });
        cube c=new cube();   
        f.add(c);
        f.setVisible(true);
        f.setResizable(false);
	}

}
