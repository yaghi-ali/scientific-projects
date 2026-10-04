"""Portable replay of IN/OUT logic on existing track IDs, not an object detector.

Public reconstruction of the counting stage of the STM32N6/ESP32 project.
No vendor firmware, learned weights, faces or site recordings are included.
"""
from dataclasses import dataclass, field
import math


@dataclass
class LineCounter:
    line_y: float = .5
    hysteresis: float = .03
    max_gap: int = 10
    entered: int = 0
    exited: int = 0
    tracks: dict = field(default_factory=dict)
    last_frame: int = -1

    def update(self, frame, positions):
        """positions maps stable track ID to normalized y. Downward = IN.

        Frames must increase. IDs absent for more than max_gap frames expire.
        Dead-band positions preserve the last confirmed side.
        """
        if frame <= self.last_frame or self.max_gap < 0 or self.hysteresis < 0:
            raise ValueError('Invalid frame sequence or settings.')
        if any(not math.isfinite(y) or not 0 <= y <= 1 for y in positions.values()):
            raise ValueError('Positions must be finite normalized coordinates.')
        self.last_frame = frame
        self.tracks = {k:v for k,v in self.tracks.items() if frame-v[1] <= self.max_gap}
        events = []
        for track, y in positions.items():
            side = -1 if y < self.line_y-self.hysteresis else 1 if y > self.line_y+self.hysteresis else 0
            previous, _ = self.tracks.get(track, (0, frame))
            if side and previous and side != previous:
                direction = 'IN' if side == 1 else 'OUT'
                self.entered += direction == 'IN'; self.exited += direction == 'OUT'
                events.append(dict(frame=frame, track_id=track, direction=direction))
            self.tracks[track] = (side or previous, frame)
        return events


if __name__ == '__main__':
    import json
    counter = LineCounter()
    for frame, y in enumerate([.2, .4, .49, .51, .6, .8, .6, .5, .4]):
        for event in counter.update(frame, {1:y}):
            print(json.dumps(event))
    print(dict(entered=counter.entered, exited=counter.exited))
