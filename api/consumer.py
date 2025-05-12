import json
from channels.generic.websocket import AsyncWebsocketConsumer
from asgiref.sync import async_to_sync

class ProgressConsumer(AsyncWebsocketConsumer):

    async def connect(self):

        self.room_name = "progress_room"
        self.room_group_name = f"progress_{self.room_name}"
        await self.channel_layer.group_add(self.room_group_name, self.channel_name)
        await self.accept()

        await self.send(text_data = json.dumps({
            'type':'connection_established',
            'message':'You are now connected'
        }))



    async def disconnect(self, close_code):
        await self.channel_layer.group_discard(self.room_group_name, self.channel_name)

    async def receive(self, text_data):
        return None

    
    async def send_progress(self, text_data):

        print(type(text_data))

        await self.send(text_data = json.dumps({
            "progress":text_data['progress_data']['progress'],
            "status":text_data['progress_data']["status"],
        }))