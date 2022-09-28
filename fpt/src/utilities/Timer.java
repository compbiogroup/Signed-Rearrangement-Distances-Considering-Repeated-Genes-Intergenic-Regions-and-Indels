package utilities;

public class Timer {
    final long MAX_RUNTIME=600;
    long startTime;
    long mark;

    public Timer () {
        startTime = System.currentTimeMillis();
        mark = startTime;
    }

    public void mark_time() {
        mark = System.currentTimeMillis();
    }

    public double since_last_mark() {
        long current = System.currentTimeMillis();
        long elapsed = current - mark;
        return (double) elapsed / 1000;
    }

    public double remaning_time() {
        return Math.max(0.0, MAX_RUNTIME - elapsed_time());
    }

    public double elapsed_time() {
        long current = System.currentTimeMillis();
        long elapsed = current - startTime;
        return (double) elapsed / 1000;
    }

    public Boolean done() {
        if (elapsed_time() >= MAX_RUNTIME) return true;
        else return false;
    }
}
